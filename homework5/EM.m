function obj=EM(D,K,mu0,sigma0,p,max_iter,min_inc)

mu=mu0;
sigma=sigma0;
N=size(D,2);


for iter=1:max_iter
    mu0=mu;
    P=zeros(K,N);
    for k=1:K
       
        P(k,:)=Dmvnpdf(D,mu(:,k),diag(sigma(:,k)));
    end
    maxP=max(P,[],1);
    P=exp(P-repmat(maxP,K,1));
    
    h=diag(p)*P*diag(1./sum(diag(p)*P,1)); %pitfall
    
    jump=0;
    for k=1:K
        if ( sum(h(k,:))<eps*100000)
%             k1=floor(rand()*(K-1)+1);
%             while(k1==k)
%                 k1=floor(rand()*(K-1)+1);
%             end
%             mu(:,k)=rand(size(D,1),1)*0.01-0.02+mu(:,k1);
%             sigma(:,k)=ones(size(D,1),1);
            
                 mu_ind=randperm(size(D,2));
                 mu_ind=mu_ind(1);
                 mu(:,k)=D(:,mu_ind)';
                 sigma(:,k)=var(D,1,2);
            p(k)=sum(h(k,:))/N;
            disp('disturb to avoid local optima');
            jump=1;
            continue;
        end
        mu(:,k)=sum(D*diag(h(k,:)),2)/sum(h(k,:));
        sigma(:,k)=sum((D-repmat(mu(:,k),1,N)).^2*diag(h(k,:)),2)/sum(h(k,:));
        sigma(:,k)=sigma(:,k)+0.0001*ones(size(D,1),1);
        if (any((sigma(:,k))<1e-7))
%             k1=floor(rand()*(K-1)+1);
%             while(k1==k)
%                 k1=floor(rand()*(K-1)+1);
%             end
%             mu(:,k)=rand(size(D,1),1)*0.01-0.02+mu(:,k1);
%             sigma(:,k)=ones(size(D,1),1);
                 mu_ind=randperm(size(D,2));
                 mu_ind=mu_ind(1);
                 mu(:,k)=D(:,mu_ind)';
                 sigma(:,k)=var(D,1,2);
            p(k)=sum(h(k,:))/N;
            disp('disturb to avoid local optima');
            jump=1;
            continue;
        end
%         
%          [~,err] = cholcov(diag(sigma(:,k)),0);
%             if err ~= 0
%                 disp(diag(sigma(:,k)));
%                 disp(h(k,:));
%                 error(message('stats:mvnpdf:BadMatrixSigma'));
%                 
%             end
        p(k)=sum(h(k,:))/N;
    end
    if (jump==1)
        p=rand(K,1);
        p=p/sum(p);
    end
    if(norm(mu-mu0,'fro')<min_inc*K*size(D,1)&&iter>30&&jump==0)
        break;
    end
end
obj.mu=mu;
obj.sigma=sigma;
obj.p=p;
obj.iter=iter;
obj.update=norm(mu-mu0,'fro')/K/size(D,1);

function P=Dmvnpdf(D,mu,sigma)
N=size(D,2);
d=size(D,1);
try
A=sigma^(-1/2);
catch 
    disp(sigma)
    error('wrong');
end
D_bar=A*(D-repmat(mu,1,N));

%P=exp(-1/2*sum(D_bar.^2)-d/2*log(2*pi)-1/2*sum(log(diag(sigma))));
P=(-1/2*sum(D_bar.^2,1)-d/2*log(2*pi)-1/2*sum(log(diag(sigma))));
