clear;clc;
allK=[1 2 4 8 16 32];
alld=[1 2 4 8 16 24 32 40 48 56 64];
for K=allK
    if exist(['gmm_BG_FG_' num2str(K) '.mat'],'file')==2
        continue;
    else
        
        load TrainingSamplesDCT_8_new.mat
        
        
        
        gmm_BG=cell(1,5);
        
        %for td=1:length(alld)
            d=64;
            for iter=1:5
                D=TrainsampleDCT_BG(:,1:d);
                %random intialization
%                 min_D=min(D);
%                 max_D=max(D);
%                 mu=diag(max_D-min_D)*randn(d,K)+repmat(min_D',1,K);
%                 sigma=diag(max_D-min_D+1)*0.5*(ones(d,K)+rand(d,K));

                %random permutation
                mu_ind=randperm(size(D,1));
                mu_ind=mu_ind(1:K);
                mu=D(randperm(mu_ind,:))';
                sigma=repmat(var(D)',1,K);
                p=rand(K,1);
                p=p/sum(p);
                gmm_BG{iter}=EM(D',K,mu,sigma,p,2000,0.001);
                disp([';iter ' num2str(iter) '/5 ' num2str(gmm_BG{iter}.iter) ' ' num2str(gmm_BG{iter}.update)]);
            end
        %end
        
        gmm_FG=cell(1,5);
        
        %for td=1:length(alld)
            d=64;
            for iter=1:5
                D=TrainsampleDCT_FG(:,1:d);
                %random intialization
%                 min_D=min(D);
%                 max_D=max(D);
%                 mu=diag(max_D-min_D)*randn(d,K)+repmat(min_D',1,K);
%                 sigma=diag(max_D-min_D+1)*0.5*(ones(d,K)+rand(d,K));
                
                %random permutation
                
                mu_ind=randperm(size(D,1));
                mu_ind=mu_ind(1:K);
                mu=D(randperm(mu_ind,:))';
                sigma=repmat(var(D)',1,K);
                p=rand(K,1);
                p=p/sum(p);
                gmm_FG{iter}=EM(D',K,mu,sigma,p,2000,0.001);
                disp([';iter ' num2str(iter) '/5 ' num2str(gmm_FG{iter}.iter) ' ' num2str(gmm_FG{iter}.update)]);
            end
        %end
        p_BG=size(TrainsampleDCT_BG,1);
        p_FG=size(TrainsampleDCT_FG,1);
        p_BG=p_BG/(p_BG+p_FG);
        p_FG=1-p_BG;
        
        save(['gmm_BG_FG_' num2str(K) '.mat'],'gmm_BG','gmm_FG','p_BG','p_FG');
        disp(['done ' num2str(K) ' components'])
    end
end
%%
load testFeature.mat

mask=imread('cheetah_mask.bmp');
alld=[1 2 4 8 16 24 32 40 48 56 64];
allK=[1 2 4 8 16 32];
errors=cell(1,length(allK));
for tK=1:length(allK)
    K=allK(tK);
    load(['gmm_BG_FG_' num2str(K) '.mat']);
    error=zeros(length(alld),5,5);
    for td=1:length(alld)
        for t_BG=1:5
            for t_FG=1:5
                [result,error(td,t_BG,t_FG)]=evalClass_GMM(testFeature(:,1:alld(td))',K,alld(td),p_BG,gmm_BG{t_BG},p_FG,gmm_FG{t_FG},mask);
                disp([num2str(tK) '/' num2str(length(allK)) ';' num2str(td) '/' num2str(length(alld)) ';iter ' num2str(t_BG*5-5+t_FG) '/25 ']);
                %imshow(result);
                %pause(0.5);
                %imwrite(result,['result_K' num2str(all(tK)) '_D' num2str(alld(td)) '_B' num2str(t_BG) '_F' num2str(t_FG) '_' num2str(error(td,t_BG,t_FG)) '.bmp']);
                
            end
        end
    end
    errors{1,tK}=error;
end
save('errors.mat','errors');

%%
alld=[1 2 4 8 16 24 32 40 48 56 64];
load('errors.mat');
plot(alld,reshape(errors{4},11,25),'LineWidth',2)
set(gca,'FontSize',13,'FontWeight','Bold');
xlim([1,64]);
title('error rate - dimension(25 classifiers)','FontSize',15,'FontWeight','Bold');
saveas(gca, 'error_8.eps','epsc');

X=zeros(11,6);
for i=1:6
X(:,i)=mean(reshape(errors{i},11,25),2);end;
plot(alld,X,'LineWidth',2)
set(gca,'FontSize',13,'FontWeight','Bold');
title('error rate - dimension(C=mixture numbers)','FontSize',15,'FontWeight','Bold');
h_l=legend('C=1','C=2','C=4','C=8','C=16','C=32');
set(h_l,'FontSize',14,'FontWeight','Bold');
xlim([1,64]);
saveas(gca, 'errors.eps','epsc');
