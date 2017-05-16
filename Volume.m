function [ v ] = Volume( structure,proportions )

% A rough estimate of the volume of the computed structures based on PCA.
% Suggests of better estimating methods are higly appreciated.

secn=10;
v=0;
for m=1:size(structure,3)

    [~,s,~]=pca(structure(:,:,m));
    ss=sortrows(s,1);
    L=ss(end,1)-ss(1,1);
    sb=1;
    sec=1;
    sv=0;
    for i=1:size(ss,1)
        if ss(i,1)-ss(1,1)>=L/secn*sec
            sec=sec+1;
            se=i;
            t=max(ss(sb:se,:))-min(ss(sb:se,:));
            sv=sv+t(3)*t(2)*L/secn;
            sb=se;
        end
    end
    
    v=v+proportions(m)*sv;
end

end

