function [ X,loci,sizepara ,index] = Preprocess(HiC_file,loci_file,input_sizepara)

% Data Preprocess before modeling.


% Data input
X=load(HiC_file);
loci=load(loci_file);
n=size(loci,1);
delt=loci(2:n)-loci(1:n-1);
[delt,sortind]=sort(delt,'descend');

% Compute the packing density (see "Methods" in paper)
% [10E3,50E3,100E3,250E3,500E3,1E6]: Bins for some common resolutions (10k,50k,100k,150k,500k,1M). 
% Users can change them according to the requirement.
mass_density=[10E3,50E3,100E3,250E3,500E3,1E6]./([10E3,50E3,100E3,250E3,500E3,1E6].^(1/3)*1E4/130/1E4^(1/3));
commonres=mode(delt);
sizepara=mass_density(ResolutionMap(commonres));

if input_sizepara > 0
    sizepara=input_sizepara;
end
% Build the "beads-on-a-string" model for the uninformed long sequences
if mass_density(ResolutionMap(delt(1)))>sizepara
    i=1;
    while  mass_density(ResolutionMap(delt(i)))>sizepara
        rec(i,1)=sortind(i);
        rec(i,2)=floor(delt(i)/commonres)-1;
        i=i+1;
    end
    rec=sortrows(rec,1);
    totaln=n+sum(rec(:,2));  % totaln: Total numbers of loci ("beads") of the "beads-on-a-string" model
    tloci=zeros(totaln,1);   % tloci: Loci of the "beads-on-a-string" model
    tX=zeros(totaln,totaln); % tX: Hi-C map of the "beads-on-a-string" model
    b=1;nc=1;index=[];       % b: Pointer for begin locus.  
                             % nc: Counter for loci of the "beads-on-a-string" model
    for i=1:size(rec,1)
        f=rec(i,1);          % f: Pointer for end locus.
        tloci(nc:nc+f-b,1)=loci(b:f,1);
        index=[index;(nc:nc+f-b)'];
        pre=loci(f);
        nc=nc+f-b+1;
        for j=1:rec(i,2)
            tloci(nc,1)=pre+commonres;
            pre=tloci(nc,1);
            nc=nc+1;
        end
        b=rec(i,1)+1;
    end
    f=n;
    tloci(nc:nc+f-b,1)=loci(b:f,1);
    index=[index;(nc:nc+f-b)'];
    
    loci=tloci;
    tX(index,index)=X;
    X=tX;
else
    index=find(loci>-1);
end

% Make sure that all diagonal numbers are zero
X=X-diag(diag(X));
end

