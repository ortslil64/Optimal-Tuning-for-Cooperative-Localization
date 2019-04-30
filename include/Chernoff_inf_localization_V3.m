function [ Ci ] = Chernoff_inf_localization_V3( PF,n,alpha,beta,Xt,K1,K2 )

PF1 = PF.X{1,n};
PF2 = PF.X{2,n};
PF3 = PF.X{3,n};
w1 = PF.W{1,n};
w2 = PF.W{2,n};
w3 = PF.W{3,n};

TF12 = Xt{1,2}(1:2);
TF13 = Xt{1,3}(1:2);
particles2rel1 = [PF2(:,1)+TF12(1),PF2(:,2)+TF12(2),PF1(:,3)];
particles3rel1 = [PF3(:,1)+TF13(1),PF3(:,2)+TF13(2),PF1(:,3)];


i1 = find(w1>0.00001);
w1 = w1(i1);
mu1 = PF1(i1,1:2);
for kk = 1:size(mu1,1)
    sigma1(:,:,kk)=beta*eye(2);
end
GM1 = gmdistribution(mu1,sigma1,w1);

% PF = [PF1;particles2rel1];
% PF1 = PF(randsample(2*Np,Np),:);


i2 = find(w2>0.00001);
w2 = w2(i2);
mu2 = particles2rel1(i2,1:2);
for kk = 1:size(mu2,1)
    sigma2(:,:,kk)=beta*eye(2);
end
GM2 = gmdistribution(mu2,sigma2,w2);

i3 = find(w3>0.00001);
w3 = w3(i3);
mu3 = particles3rel1(i3,1:2);
for kk = 1:size(mu3,1)
    sigma3(:,:,kk)=beta*eye(2);
end
GM3 = gmdistribution(mu3,sigma3,w3);

M = mean([mu1;mu2;mu3],1);
C = cov([mu1;mu2;mu3]);
GM = gmdistribution(M,C);
GM_IS = fitgmdist([mu1;mu2;mu3],K1,'RegularizationValue',0.01,'Replicates',1,'Options',statset('Display','off','MaxIter',50,'TolFun',1e-5));
warning('off','stats:gmdistribution:FailedToConvergeReps');
g_samples = random(GM_IS,K2);
w123 = ((pdf(GM1,g_samples).^(alpha(1))).*(pdf(GM2,g_samples).^(alpha(2))).*(pdf(GM3,g_samples).^(alpha(3))))./pdf(GM_IS,g_samples);
Ci = sum(w123)/K2;
end

