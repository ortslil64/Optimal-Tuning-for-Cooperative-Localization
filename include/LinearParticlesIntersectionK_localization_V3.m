function [ particles ] = LinearParticlesIntersectionK_localization_V3( PF_array,n,alpha,beta,Xt,to_i )

Np = length(PF_array.W{1,n});

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
a1 = ceil(Np*alpha(1));
a2 = ceil(Np*alpha(2));
a3 = Np - a1 - a2;
particles = [PF1(randsample(Np,a1,true,w1),:);...
    particles2rel1(randsample(Np,a2,true,w2),:);...
    particles3rel1(randsample(Np,a3,true,w3),:)];

end

