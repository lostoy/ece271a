function [result,error]=evalClass_GMM(D,K,d,p_BG,gmm_BG,p_FG,gmm_FG,mask)
[h,w]=size(mask);
N=size(D,2);
P_BG=zeros(K,N);
for k=1:K
    P_BG(k,:)=Dmvnpdf(D,gmm_BG.mu(1:d,k),diag(gmm_BG.sigma(1:d,k)));
end
P_1=gmm_BG.p'*P_BG;

P_FG=zeros(K,N);
for k=1:K
    P_FG(k,:)=Dmvnpdf(D,gmm_FG.mu(1:d,k),diag(gmm_FG.sigma(1:d,k)));
end
P_2=gmm_FG.p'*P_FG;

result=log(P_1)+log(p_BG)<log(P_2)+log(p_FG);

result=reshape(result(1:(w-7)*(h-7)),w-7,h-7)';

mask=mask(4:h-4,4:w-4);


error=sum(sum((abs(result*1.0-double(mask)/255))))/h/w;

function P=Dmvnpdf(D,mu,sigma)
N=size(D,2);
d=size(D,1);
A=sigma^(-1/2);

D_bar=A*(D-repmat(mu,1,N));

P=exp(-1/2*sum(D_bar.^2,1))/((2*pi)^(d/2)*det(sigma)^(1/2));
