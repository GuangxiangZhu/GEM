function [ ] = BayesParaSelect( begin_para, end_para, real_volume, HiC_file, loci_file, max_iter, M, infer_latent )

% Parameter selection based on Bayesian approach

% begin_para and end_para: Select parameters in the range of [5*10^begin_para, 5*10^end_para]. 
% real_volume: Real volume of the chromatin. If you do not have the priori 
%              information of the real volume, you can set real_volume -1
%              or -2 to use the estimated value provided by GEM (-1 for 
%              human cell, -2 for yeast cell).
% HiC_file: File name of Hi-C map.
% loci_file: File name of genomic loci.
% max_iter: Maximum number of iterations. 
% M: Number of conformations.
% lambdaE: Energy coefficient.
% infer_latent: Whether to infer the latent function (1/0).


% v0: The real volume from user input or estimated by GEM
loci=load(loci_file);
if real_volume==-1
    % DNA density for human nuclei : 0.012 bp/nm3
    v0=(loci(end)-loci(1))/12E6*1E9;
elseif real_volume==-2
    % DNA density for yeast nuclei : 0.006 bp/nm3
    v0=(loci(end)-loci(1))/6E6*1E9;
else
    v0=real_volume;
end

score_best=0;
for para=begin_para:end_para
    [structure_tmp,proportions_tmp,C_tmp,C1_tmp,C2_tmp]=GEMpara(HiC_file,loci_file,max_iter,M,5*10^para);
    socre_tmp=BayesianScore( structure_tmp,proportions_tmp,C1_tmp,v0 );
    if socre_tmp > score_best
        score_best=socre_tmp;
        lambdaE_best=5*10^para;
        structure_best=structure_tmp;
        proportions_best=proportions_tmp;
        C_best=C_tmp;
        C1_best=C1_tmp;
        C2_best=C2_tmp;
    end
end    

disp('============================ RESULTS ============================');
X=load(HiC_file);
if infer_latent==1
    LatentFunction( X,structure_best,M,proportions_best );
end

%The conformations all have been aligned using the singular value decomposition (SVD) algorithm
if M>1
    for i=1:M
        [ RMSE,structure_best(:,:,i) ] = SVD3D( structure_best(:,:,i),structure_best(:,:,1));
    end
end
% Output reconstructed structure to 'conformation[1-M].txt'
for m=1:M
    dlmwrite(['conformation' num2str(m) '.txt'],structure_best(:,:,m),'delimiter', '\t','precision','%6.4f');
end
% Output the corresponding weights of conformations to 'proportions.txt'
dlmwrite('proportions.txt',proportions_best,'delimiter', '\t','precision','%6.4f');

disp (['Optimal parameter (lambdaE): ' num2str(lambdaE_best)]);
disp (['Total cost (C): ' num2str(C_best)]);
disp (['Data cost (C1): ' num2str(C1_best)]);
disp (['Energy cost (C2): ' num2str(C2_best)]);
disp (['The reconstructed conformations and the conformation weights are written to ''conformation[1-' num2str(M) '].txt'' and ''proportions.txt'', respectively. ' ]);

% ------------------- Latent function  -------------------- %
% Capture the latent function by comparing the modeled      %
% structures with the original Hi-C data.                   %
% --------------------------------------------------------- %


disp('=================================================================');
end

