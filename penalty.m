function [ydata_min,cost_min,proportions_min,cost1_min,costE_min,epochcost]=penalty(P,loci,max_iter,sizepara,M,lambdaE)

% This function is for penalty method
% Optimization method is similar with the one in t-SNE(Laurens van der Maaten, 2008)
% Exterior point penalty function method for increasingly big delta to
% constrait sequence distance between two neighbored loci.
n=size(loci,1);
for runtime=1:3
    ydata=5000*rand(n,3,1);
    %ydata=one_by_one(loci,sizepara);
    ydata = bsxfun(@minus, ydata, mean(ydata, 1));
    %ydata=MDS(loci,P);
    %one conformation
    [ydata,cost,proportions,cost1,costE,epochcost1] = gradient( ydata, P,loci,max_iter,sizepara,1,lambdaE );
    %%{
    %disp(['single ' num2str(cost1) ' ' num2str(costE)]);
    ydata=ydata(:,:,ones(1,M))+noise(M)*randn(n,3,M);
    ydata = bsxfun(@minus, ydata, mean(ydata, 1));
   
    [ydata,cost,proportions,cost1,costE,epochcost2] = gradient( ydata, P,loci,max_iter,sizepara,M,lambdaE );
    disp(['multi ' num2str(cost1) ' ' num2str(costE)]);
    % For find the smallest cost and its ydata.
    if runtime==1
        cost_min=cost;
        ydata_min=ydata;
        proportions_min=proportions;
        cost1_min=cost1;
        costE_min=costE;
        epochcost=[epochcost1];
    end
    if runtime >1
        if cost_min>cost
            cost_min=cost;
            ydata_min=ydata;
            proportions_min=proportions;
            cost1_min=cost1;
            costE_min=costE;
            epochcost=[epochcost1;0;epochcost2];
        end
    end
    %}
end
end