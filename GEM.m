function []=GEM(HiC_file,loci_file,max_iter,M,lambdaE,infer_latent,input_sizepara)

% ----------------------- Overview ------------------------ %
% Main function for GEM                                     %
% --------------------------------------------------------- %
% HIC_file: File name of Hi-C map                           %
% loci_file: File name of genomic loci                      %
% max_iter: Maximum number of iterations                    %
% M: Number of conformations                                %
% lambdaE: Energy coefficient                               %
% infer_latent: Whether to infer the latent function (1/0)  %
% input_sizepara: Packing density provided by user.         %
%                 -1 means using the estimated value.       %
% --------------------------------------------------------- %


% ---------------------- Preproccess ---------------------- %
% Data Preprocess before modeling.                          %
% --------------------------------------------------------- %
% X: Hi-C map                                               %
% loci: Genomic loci                                        %
% sizepara: Packing density                                 %
% index: The index of the genomic loci of interest          %
% --------------------------------------------------------- %
[ X,loci,sizepara,index ] = Preprocess(HiC_file,loci_file,input_sizepara);


% ---------------------- Optimizer ------------------------ %
% An optimization process that considers both Hi-C data     %
% and conformation energy.                                  %
% --------------------------------------------------------- %
% structure: The reconstructed ensemble of conformations    %
%            (N*3*M matrix)                                 %
% proportions: The corresponding weights of conformations   %
% C: Total cost (C1+lambda_E*C2)                            %
% C1: Data cost (KL divergence)                             %
% C2: Energy cost (Conformation energy)                     %
% --------------------------------------------------------- %
[structure,proportions,C,C1,C2]=Optimizer(X,loci,max_iter,sizepara,M,lambdaE,index);


disp('============================ RESULTS ============================');

%The conformations all have been aligned using the singular value decomposition (SVD) algorithm
if M>1
    for i=1:M
        [ RMSE,structure(:,:,i) ] = SVD3D( structure(:,:,i),structure(:,:,1) );
    end
end
% Output reconstructed structure to 'conformation[1-M].txt'
for m=1:M
    dlmwrite(['conformation' num2str(m) '.txt'],structure(:,:,m),'delimiter', '\t','precision','%6.4f');
end

% Output the corresponding weights of conformations to 'proportions.txt'
dlmwrite('proportions.txt',proportions,'delimiter', '\t','precision','%6.4f');

disp (['Total cost (C): ' num2str(C)]);
disp (['Data cost (C1): ' num2str(C1)]);
disp (['Energy cost (C2): ' num2str(C2)]);
disp (['The reconstructed conformations and the conformation weights are written to ''conformation[1-' num2str(M) '].txt'' and ''proportions.txt'', respectively. ' ]);

% ------------------- Latent function  -------------------- %
% Capture the latent function by comparing the modeled      %
% structures with the original Hi-C data.                   %
% --------------------------------------------------------- %
if infer_latent==1
    LatentFunction( X(index,index),structure,M,proportions );
end

disp('=================================================================');



end
