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
                min_D=min(D);
                max_D=max(D);
                mu=diag(max_D-min_D)*randn(d,K)+repmat(min_D',1,K);
                sigma=diag(max_D-min_D+1)*0.5*(ones(d,K)+rand(d,K));
                p=rand(K,1);
                p=p/sum(p);
                option=statset('MaxIter',2000);
                gmm_BG{iter}=gmdistribution.fit(D,K,'CovType','diagonal','Regularize',0.001,'Options',option);
                disp([';iter ' num2str(iter) '/5 ']);
            end
        %end
        
        gmm_FG=cell(1,5);
        
        %for td=1:length(alld)
            d=64;
            for iter=1:5
                D=TrainsampleDCT_FG(:,1:d);
                %random intialization
                min_D=min(D);
                max_D=max(D);
                mu=diag(max_D-min_D)*randn(d,K)+repmat(min_D',1,K);
                sigma=diag(max_D-min_D+1)*0.5*(ones(d,K)+rand(d,K));
                p=rand(K,1);
                p=p/sum(p);
                option=statset('MaxIter',2000);
                gmm_FG{iter}=gmdistribution.fit(D,K,'CovType','diagonal','Regularize',0.001,'Options',option);
                disp([';iter ' num2str(iter) '/5 ' ]);
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
                [~,error(td,t_BG,t_FG)]=evalClass_GMM_gmm(testFeature(:,1:alld(td))',K,alld(td),p_BG,gmm_BG{t_BG},p_FG,gmm_FG{t_FG},mask);
                disp([num2str(td) '/' num2str(length(alld)) ';iter ' num2str(t_BG*5-5+t_FG) '/25 ']);
                %imshow(result);
                %pause(0.5);
                
                
            end
        end
    end
    errors{1,tK}=error;
end
X=zeros(11,6);
for i=1:6
X(:,i)=mean(reshape(errors{i},11,25),2);end;
plot(X)
saveas(gca, 'errors.eps','epsc');
save('errors.mat','errors');
save('errors.mat','errors');