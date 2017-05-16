function [ydata,proportions,C,C1,C2] = BasicOptCell(ydata, P,loci,max_iter,sizepara,M,lambdaE   )

% Basic optimizer cell
% If M = 1, it serves as an optimizer for average structure
% If M > 1, it serves as an optimizer for multiple conformations
% Some gradients refer to t-SNE (Laurens van der Maaten, 2008)
% See "Methods" in the paper of GEM for more 
% Suggests of better optimization are higly appreciated.

% Joint probability P in Hi-C space
P=P/sum(sum(P));
zeroP=find(P==0);         % Zero interactions
unzeroP=find(P>0);        % Non-zero iteractions
mmax=1E5;                 % Prevent the excessive value  

% Coeffiencts
epsilonY = 100;           % Learning rate for Y
epsilonW =0;              % Learning rate for W (initialization)
changeflag=0;             % Whether epsilonW is changed
minerr=10^(-6);           % Optimization stops when cost function is less than minerr for some time which is counted by qcount                           
qcount=0;                                 

% Compute the spatial distance between two continuous loci
n=size(loci,1);
l=loci(2:n)-loci(1:n-1);
l=l/sizepara ;

% Initialization
weights=zeros(M,1);
proportions = exp(-weights);
proportions= proportions/sum(proportions);
dCdP = zeros(M,1);
dCdD = zeros(n,n,M);
dCdY = zeros(size(ydata));

% Optimization begins
for iter=1:max_iter
    
    % Compute the conformation energy and its gradient
    for m=1:M
        [ E,gE ] = Energy( ydata(:,:,m),l,n);
        tmpE(m)=E;                     
        tmpgE(:,:,m)=gE;
    end
    
    % Compute pairwise affinities per conformation
    for m=1:M
    	sum_ydata = sum(ydata(:,:,m) .^ 2, 2);
    	tmp = 1 ./ (1 + bsxfun(@plus, sum_ydata, bsxfun(@plus, sum_ydata', -2 * ydata(:,:,m) * ydata(:,:,m)')));
    	tmp(zeroP)=0;
        tmp(1:n+1:end) = 0;                          % set diagonal to zero
        num(:,:,m) = tmp;
    end
        
    % Compute joint probability that point i and j are neighbors
    QZ=zeros(n,n);
    for m=1:M
        QZ = QZ + proportions(m)*num(:,:,m);
    end
    Z = sum(QZ(:));        
    Q = QZ ./ Z;
    
    % Compute the gradient of C1 with respect to mixture proportions
    QP = Q - P;
    for m=1:M
        tmp=(1 ./ QZ) .* QP.*num(:,:,m);
        dCdP(m) = sum(sum(tmp(unzeroP)))+lambdaE*tmpE(m);
    end

    % Compute the gradient of C1 with respect to mixture weights
    dCdW = proportions .* ( sum(dCdP .* proportions)-dCdP);

    % Compute the gradient of C1 with respect to pairwise distances
    for m=1:M
    	tmp = proportions(m)*( -QP) .* num(:,:,m).* num(:,:,m)./ QZ;
        tmp(zeroP)=0;
        dCdD(:,:,m)=tmp;
    end
    
    % Compute the gradient of C with respect to the coordinate
    for m=1:M
        for i=1:n
            dCdY(i,:,m) = 4*sum(bsxfun(@times, dCdD(:,i,m),bsxfun(@minus, ydata(i,:,m), ydata(:,:,m))), 1);
        end
        dCdY(:,:,m)=dCdY(:,:,m)+lambdaE*proportions(m)*tmpgE(:,:,m);
    end
    
    % Cost
    % C: Total cost (C1+lambda_E*C2)                           
    % C1: Data cost (KL divergence)                           
    % C2: Energy cost (Conformation energy)  
    C2=0;
    for m=1:M
        C2=C2+proportions(m)*tmpE(m);
    end 
    C = sum(P(unzeroP) .* log(P(unzeroP) ./ Q(unzeroP)))+lambdaE*C2;
    C1=C-lambdaE*C2;
    if iter==1
        pC=C;
    end
    
    % Find a local optimal solution using adaptive learning rate
    if  C>pC
        epsilonY=epsilonY*0.9;
        epsilonW=epsilonW*0.9;
        ydata=pydata;
        weights=pweights;
        proportions=pproportions;
        continue;
    else
        if pC-C<minerr
            qcount=qcount+1;
        else
            qcount=0;
        end
        
        if~rem(iter, 100)
            epsilonY=epsilonY*2;
        end
        pC=C;
    end
    
    if qcount>10 && epsilonY<1
        break;
    end
    
    if ~rem(iter, 100)
    %    disp(['Iteration ' num2str(iter) ': C is ' num2str(C) ', C1 is ' num2str(C1) ', C2 is ' num2str(C2) ]);
    end
    
    if C1<0.2
        changeflag=1;
    end
    if changeflag
        epsilonW=1;
    end
    
    pydata=ydata;
    pweights=weights;
    pproportions=proportions;
    if iter~=max_iter    
        ydata=ydata-epsilonY*dCdY;
        ydata = bsxfun(@minus, ydata, mean(ydata, 1));
        weights = weights - epsilonW * dCdW;
        proportions = exp(-weights);
        proportions(proportions>mmax)=mmax;
        proportions= proportions/sum(proportions);
    end
end
end