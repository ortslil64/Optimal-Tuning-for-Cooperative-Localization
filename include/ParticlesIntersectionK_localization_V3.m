function [ w123 ] = ParticlesIntersectionK_localization_V3( PF_array,n,alpha,beta,Xt,to_i )



if to_i==1
    w1 = PF_array.W{1,n};
    w2 = PF_array.W{2,n};
    w3 = PF_array.W{3,n};
    PF1 = PF_array.X{1,n};
    PF2 = PF_array.X{2,n};
    PF3 = PF_array.X{3,n};
    TF12 = Xt{1,2}(1:2);
    TF13 = Xt{1,3}(1:2);
    particles2rel1 = [PF2(:,1)+TF12(1),PF2(:,2)+TF12(2),PF1(:,3)];
    particles3rel1 = [PF3(:,1)+TF13(1),PF3(:,2)+TF13(2),PF1(:,3)];
end
if to_i==2
    w1 = PF_array.W{2,n};
    w2 = PF_array.W{1,n};
    w3 = PF_array.W{3,n};
    PF1 = PF_array.X{2,n};
    PF2 = PF_array.X{1,n};
    PF3 = PF_array.X{3,n};
    TF12 = Xt{2,1}(1:2);
    TF13 = Xt{2,3}(1:2);
    particles2rel1 = [PF2(:,1)+TF12(1),PF2(:,2)+TF12(2),PF1(:,3)];
    particles3rel1 = [PF3(:,1)+TF13(1),PF3(:,2)+TF13(2),PF1(:,3)];
end
if to_i==3
    w1 = PF_array.W{3,n};
    w2 = PF_array.W{2,n};
    w3 = PF_array.W{1,n};
    PF1 = PF_array.X{3,n};
    PF2 = PF_array.X{2,n};
    PF3 = PF_array.X{1,n};
    TF12 = Xt{3,2}(1:2);
    TF13 = Xt{3,1}(1:2);
    particles2rel1 = [PF2(:,1)+TF12(1),PF2(:,2)+TF12(2),PF1(:,3)];
    particles3rel1 = [PF3(:,1)+TF13(1),PF3(:,2)+TF13(2),PF1(:,3)];
end
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


w123 = (pdf(GM1,PF1(:,1:2)).^alpha(1)).*(pdf(GM2,PF1(:,1:2)).^alpha(2)).*(pdf(GM3,PF1(:,1:2)).^alpha(3));
end

