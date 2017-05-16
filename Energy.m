function [ E,gE ] = Energy( s,l,n)



% Parameters for energy determination
ks=2.3*10^-21;
ktheta=1.7*10^-20;
KT=4.1*10^-21;
delta=30;
tmatrix=triu(ones(n,n),1);
tmatrix=tmatrix([2:end 1],:);
tmatrix(end,:)=0;

% Conformation Energy
% Es: Stretching energy
s1=s(1:n-1,:);
s2=s(2:n,:);
ds=s2-s1;
normds=sqrt(sum((ds).^2,2));
errds=normds-l;
Es=ks/2*sum(errds.^2);

% Eb: Bending energy
ds1=ds(1:n-2,:);
ds2=ds(2:n-1,:);
normds1=sqrt(sum((ds1).^2,2));
normds2=sqrt(sum((ds2).^2,2));
tds2=bsxfun(@rdivide,bsxfun(@rdivide,ds2,normds1),normds2);
Eb=ktheta/2*sum(acos(sum(ds1.*tds2,2)).^2);

% Ee: Excluding energy
dm=squareform(pdist(s));
dmat=dm.*tmatrix;
flag=sum(sum(dmat<1.1225*delta & dmat>0));
if flag==0
    Ee=0;
else
    tdmat=zeros(size(dmat));
    tdmat(dmat<1.1225*delta & dmat>0)=delta./dmat(dmat<1.1225*delta & dmat>0);
    Ee=4*KT*(sum(sum((tdmat.^12-tdmat.^6)))+flag/4);
end

% Total Energy
E=Es+Eb+Ee;

% Gradient of conformation energy
% gEs: Gradient of stretching energy
gEs(1,:)=-ks*( ds(1,:)*errds(1)/normds(1));
    
gEs(2:n-1,:)=ks*( -bsxfun(@times,bsxfun(@rdivide,ds(2:n-1,:),normds(2:n-1)),errds(2:n-1))...
        + bsxfun(@times,bsxfun(@rdivide,ds(1:n-2,:),normds(1:n-2)),errds(1:n-2)) );

gEs(n,:)=ks*( ds(n-1,:)*errds(n-1)/normds(n-1) );

% gEb: Gradient of bending energy
t1=sum(ds1.*ds2,2);
t2=bsxfun(@rdivide,ds,normds.^2);
temp=(normds1.*normds2).^2-t1.^2;
temp(temp<0)=0;
t3=sqrt(temp);

t3(t3==0)=10^-10;

temp=t1./normds1./normds2;
temp(temp<0)=0;
temp(temp>1)=1;
t4=ktheta*acos(temp)./t3;
G1=zeros(size(s));
G2=zeros(size(s));
G3=zeros(size(s));
G1(2:n-1,:)=bsxfun(@times,(ds2-ds1+bsxfun(@times,t2(2:n-1,:),t1)-bsxfun(@times,t2(1:n-2,:),t1)),t4);
G2(3:n,:)=bsxfun(@times,(ds1-bsxfun(@times,t2(2:n-1,:),t1)),t4);
G3(1:n-2,:)=bsxfun(@times,(-ds2+bsxfun(@times,t2(1:n-2,:),t1)),t4);
gEb=G1+G2+G3;


% gEe: Gradient of excluding energy
flagmat=dmat<1.1225*delta & dmat>0;
for i=1:n
    gEe(i,:)=[0 0 0];
    for j=1:n
        if flagmat(i,j)~=0
            gEe(i,:)=gEe(i,:)+(-2*delta^12/dm(i,j)^14+delta^6/dm(i,j)^8)*(s(i,:)-s(j,:));
        end
    end
end
gEe=24*KT*gEe;

% Total gradient of conformation energy
gE=gEs+gEb+gEe;
end

