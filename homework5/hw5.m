D=[-2.5 -2 -1 2.5 2;-1 0.5 0 -1 1];
sigma=[0.1 1 2];
miu=-5:0.1:5;

loglikelihood=zeros(length(miu),length(miu));
for s=sigma
    for t1=1:length(miu)
        for t2=1:length(miu)
            m1=miu(t1);
            m2=miu(t2);
            centered_D_1=D-repmat([m1;m2],1,size(D,2));
            P_1=exp(-sum(centered_D_1.^2)/2/s)/(2*pi*s);
            
            centered_D_2=D+repmat([m1;m2],1,size(D,2));
            P_2=exp(-sum(centered_D_2.^2)/2/s)/(2*pi*s);
            
            loglikelihood(t1,t2)=sum(log(P_1/2+P_2/2));
            
        end
    end
    surf(miu,miu,loglikelihood);
    set(gca,'FontSize',13,'FontWeight','Bold');
    title(['loglikelihood, sigma=' num2str(s)],'FontSize',15,'FontWeight','Bold');
    saveas(gca, ['loglikelihood_' num2str(s) '.eps'],'epsc');
    pause;
    close all;
end

%% EM
D=[-2.5 -2 -1 2.5 2;-1 0.5 0 -1 1];

p_1=0.5;
p_2=0.5;
sigma1_1=1;
sigma1_2=1;
sigma2_1=1;
sigma2_2=1;
mu1=[-.1,0]';
mu2=-mu1;
M=300;
mu1s=zeros(2,M+1);
mu2s=zeros(2,M+1);
mu1s(:,1)=mu1;
mu2s(:,1)=mu2;
for step=2:M+1
    
    centered_D_1=D-repmat(mu1,1,size(D,2));
    centered_D_1(1,:)=centered_D_1(1,:)/sqrt(sigma1_1);
    centered_D_1(2,:)=centered_D_1(2,:)/sqrt(sigma1_2);
    P_1=exp(-sum(centered_D_1.^2)/2)/(2*pi*sqrt(sigma1_1*sigma1_2));
    
    centered_D_2=D-repmat(mu2,1,size(D,2));
    centered_D_2(1,:)=centered_D_2(1,:)/sqrt(sigma2_1);
    centered_D_2(2,:)=centered_D_2(2,:)/sqrt(sigma2_2);
    P_2=exp(-sum(centered_D_2.^2)/2)/(2*pi*sqrt(sigma2_1*sigma2_2));
    
    h1=P_1*p_1./(P_1*p_1+P_2*p_2);
    
    [X,Y]=meshgrid(-4:0.1:4);
    X=[X(:)';Y(:)'];
    
    centered_X_1=X-repmat(mu1,1,size(X,2));
    centered_X_1(1,:)=centered_X_1(1,:)/sqrt(sigma1_1);
    centered_X_1(2,:)=centered_X_1(2,:)/sqrt(sigma1_2);
    P_1=exp(-sum(centered_X_1.^2)/2)/(2*pi*sqrt(sigma1_1*sigma1_2));
    
    centered_X_2=X-repmat(mu2,1,size(X,2));
    centered_X_2(1,:)=centered_X_2(1,:)/sqrt(sigma2_1);
    centered_X_2(2,:)=centered_X_2(2,:)/sqrt(sigma2_2);
    P_2=exp(-sum(centered_X_2.^2)/2)/(2*pi*sqrt(sigma2_1*sigma2_2));
    hX_1=P_1*p_1./(P_1*p_1+P_2*p_2);
    
    
    
    mu1=sum(D*diag(h1),2)/sum(h1);
    mu2=sum(D*diag(1-h1),2)/sum(1-h1);
    
    sigma1_1=sum((D(1,:)-mu1(1)).^2*diag(h1))/sum(h1);
    sigma1_2=sum((D(2,:)-mu1(2)).^2*diag(h1))/sum(h1);
    
    sigma2_1=sum((D(1,:)-mu2(1)).^2*diag(1-h1))/sum(1-h1);
    sigma2_2=sum((D(2,:)-mu2(2)).^2*diag(1-h1))/sum(1-h1);
    
    p_1=sum(h1)/size(D,2);
    p_2=1-p_1;
    
    mu1s(:,step)=mu1;
    mu2s(:,step)=mu2;
    
    if (step<=4 )
    
     figure(1)
     plot(mu1s(1,1:step),mu1s(2,1:step),'r+',mu2s(1,1:step),mu2s(2,1:step),'b*','MarkerSize',12);
     xlim([-3,3]);
     ylim([-3,3]);
     hold on;
     plot(D(1,:),D(2,:),'gs','MarkerSize',12);
     
     el_1_x=-sqrt(sigma1_1):0.01:sqrt(sigma1_1);
     el_1_y=sqrt(1-el_1_x.^2/sigma1_1)*sqrt(sigma1_2);
     el_1_x=[el_1_x fliplr(el_1_x)]+mu1(1);
     el_1_y=[el_1_y fliplr(-el_1_y)]+mu1(2);
     
     el_2_x=-sqrt(sigma2_1):0.01:sqrt(sigma2_1);
     el_2_y=sqrt(1-el_2_x.^2/sigma2_1)*sqrt(sigma2_2);
     el_2_x=[el_2_x fliplr(el_2_x)]+mu2(1);
     el_2_y=[el_2_y fliplr(-el_2_y)]+mu2(2);
     
     plot(el_1_x,el_1_y,'r',el_2_x,el_2_y,'b:','LineWidth',2);
     
     h_l=legend('mean_1','mean_2','data','el_1','el_2');
     set(h_l,'FontSize',14,'FontWeight','Bold');
     set(gca,'FontSize',13,'FontWeight','Bold');
     title(['EM iter: ' num2str(step-1)],'FontSize',15,'FontWeight','Bold');
     
     saveas(gca, ['EM_' num2str(step-1) '.eps'],'epsc');
     
     figure(2)
     surf(-4:0.1:4,-4:0.1:4,reshape(hX_1,sqrt(length(hX_1)),sqrt(length(hX_1))));
     set(gca,'FontSize',13,'FontWeight','Bold');
     title(['p(z=1|x) EM iter: ' num2str(step-1)],'FontSize',15,'FontWeight','Bold');
     
     saveas(gca, ['pos_EM_' num2str(step-1) '.eps'],'epsc');

     pause();
     close all;
     
    end
    
    if (norm(mu1s(:,step)-mu1s(:,step-1))<=0.001 && norm(mu2s(:,step)-mu2s(:,step-1))<=0.001)
        break;
    end
end 


    figure(1)
     plot(mu1s(1,1:step),mu1s(2,1:step),'r+',mu2s(1,1:step),mu2s(2,1:step),'b*','MarkerSize',12);
     xlim([-3,3]);
     ylim([-3,3]);
     hold on;
     plot(D(1,:),D(2,:),'gs','MarkerSize',12);
     
     el_1_x=-sqrt(sigma1_1):0.01:sqrt(sigma1_1);
     el_1_y=sqrt(1-el_1_x.^2/sigma1_1)*sqrt(sigma1_2);
     el_1_x=[el_1_x fliplr(el_1_x)]+mu1(1);
     el_1_y=[el_1_y fliplr(-el_1_y)]+mu1(2);
     
     el_2_x=-sqrt(sigma2_1):0.01:sqrt(sigma2_1);
     el_2_y=sqrt(1-el_2_x.^2/sigma2_1)*sqrt(sigma2_2);
     el_2_x=[el_2_x fliplr(el_2_x)]+mu2(1);
     el_2_y=[el_2_y fliplr(-el_2_y)]+mu2(2);
     
     plot(el_1_x,el_1_y,'r',el_2_x,el_2_y,'b:','LineWidth',2);
     
     h_l=legend('mean_1','mean_2','data','el_1','el_2');
     set(h_l,'FontSize',14,'FontWeight','Bold');
     set(gca,'FontSize',13,'FontWeight','Bold');
     title(['EM iter: ' num2str(step-1)],'FontSize',15,'FontWeight','Bold');
     saveas(gca, ['EM_' num2str(step-1) '.eps'],'epsc');
     
     figure(2)
     surf(-4:0.1:4,-4:0.1:4,reshape(hX_1,sqrt(length(hX_1)),sqrt(length(hX_1))));
     set(gca,'FontSize',13,'FontWeight','Bold');
     title(['p(z=1|x) EM iter: ' num2str(step-1)],'FontSize',15,'FontWeight','Bold');
     saveas(gca, ['pos_EM_' num2str(step-1) '.eps'],'epsc');
     
     pause();
     close all;
     
%% parzen-window

a=1;
x=-a:0.01:2*a;
y=zeros(1,length(x));

h=1;

for t=1:length(x)
    tx=x(t);
    if (tx<0)
        y(t)=0;
    else
        if (tx>=0&&tx<=a)
            y(t)=1/a*(1-exp(-tx/h));
        else
            y(t)=1/a*(exp(a/h)-1)*exp(-tx/h);
        end
    end
    
end
plot(x,y,'LineWidth',2);
xlim([-a,2*a]);
ylim([-.2,1.2]);
set(gca,'FontSize',13,'FontWeight','Bold');
title('$$\bar{p}(x), h=1$$','FontSize',15,'FontWeight','Bold','Interpreter','Latex');

%%
a=1;
x=0:0.01:0.05;
y=zeros(1,length(x));

h=0.01/log(100);

for t=1:length(x)
    tx=x(t);
    if (tx<0)
        y(t)=0;
    else
        if (tx>=0&&tx<=a)
            y(t)=1/a*(1-exp(-tx/h));
        else
            y(t)=1/a*(exp(a/h)-1)*exp(-tx/h);
        end
    end
    
end
plot(x,y,'LineWidth',2);
set(gca,'FontSize',13,'FontWeight','Bold');
title('$$\bar{p}(x), h=0.0022$$','FontSize',15,'FontWeight','Bold','Interpreter','Latex');
     