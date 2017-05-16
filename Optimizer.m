function [structure_min,proportions_min,C_min,C1_min,C2_min]=Optimizer(X,loci,max_iter,sizepara,M,lambdaE,ind)

% This function is for the two-stage optimization process that considers both Hi-C data and conformation energy.

n=size(loci,1);
for trial=1:3
    disp (['Trial ' num2str(trial) ':']);
    
    disp ('Average-structure optimization');
    % Initialization
    structure=5000*rand(n,3,1);
    structure = bsxfun(@minus, structure, mean(structure, 1));
    % Average-structure optimization
    [structure,proportions,C,C1,C2] = BasicOptCell( structure, X,loci,max_iter,sizepara,1,lambdaE );
    
    disp ('Multi-conformation optimization');
    % Initialization from the precomputed average structure
    structure=structure(:,:,ones(1,M))+randn(n,3,M);
    structure = bsxfun(@minus, structure, mean(structure, 1));  
    % Multi-conformation optimization
    [structure,proportions,C,C1,C2] = BasicOptCell( structure, X,loci,max_iter,sizepara,M,lambdaE );
    
    % Find the smallest cost and its structure among three trails
    if trial==1
        C_min=C;
        structure_min=structure;
        proportions_min=proportions;
        C1_min=C1;
        C2_min=C2;
    end
    if trial >1
        if C_min>C
            C_min=C;
            structure_min=structure;
            proportions_min=proportions;
            C1_min=C1;
            C2_min=C2;
        end
    end
end

% Filter out the genomic loci of interest 
structure_min=structure_min(ind,:,:);
end