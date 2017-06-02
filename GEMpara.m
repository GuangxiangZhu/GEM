function [structure,proportions,C,C1,C2]=GEM(HiC_file,loci_file,max_iter,M,lambdaE,input_sizepara)

% ----------------------- Overview ------------------------ %
% Implement of GEM for parameter selection                  %
% --------------------------------------------------------- %
% HIC_file: File name of Hi-C map                           %
% loci_file: File name of genomic loci                      %
% max_iter: Maximum number of iterations                    %
% M: Number of conformations                                %
% lambdaE: Energy coefficient                               %
% input_sizepara: Packing density provided by user.         %
%               -1 means using the estimated value by GEM.  %
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

end
