%% clear
clc;
clear;
close all;
warning('query','all');
%% ---- Load map ---- %%
convert_map;
map_image = imread('map/map.pgm');
map_imageCropped = map_image(750:1250,950:1400);
map_imageBW = map_imageCropped < 100;
map = robotics.OccupancyGrid(map_imageBW,20);
map.GridLocationInWorld = 0.5*[-25 -25];
inflate(map,0.01);
%map.show();

clear map_image map_imageCropped map_imageBW x_size y_size vars threashold resulution rawData1 newData1 name i ii jj map_raw;
%show(map);

%% ---- Load simulation data ---- %%
load('simulation/simulation_data.mat');

%% ---- Networks and estimators parameters ---- %%
Nf = 3; % Number of robots

sigma_L = 0.5; % Likelihood param
sigma_U = 0.5; % UCLT fusion param
Np = 200; % Number of particles
Nm = 9; % Number of networks
sigma_x_step = 0.5; % process noise
sigma_y_step = 0.5; % process noise
sigma_thets_step = 0.2; % process noise
Cnoise = [0.01 0 0;0 0.01 0; 0 0 0.001]; % process noise
fusion_reg = 0.01; % fusion reglarization param
max_laser_samples = 20;
total_time_steps = 100;
start_t = 200;
MonteCarloRuns = 1; % number of Monte Carlo simulations

%% optimization
options = optimoptions('fmincon','Algorithm','interior-point','MaxIterations',100,'StepTolerance',1e-6,'Display','off');
%% initial particles
Xt = cell(Nf);
PF.X = cell(Nf,Nm);
PF.W = cell(Nf,Nm);
for ii = 1:Nf
    for jj = 1:Nm
        PF.X{ii,jj} = zeros(Np,3);
        PF.W{ii,jj} = ones(Np,1);
    end
end



for Nmc = 1:MonteCarloRuns
    start_t = randi([1 499],1,1);
    mc = mcmix(Nf);
    A=mc.P;
    [V,~] = eig(A);
    V=rand(Nf);
    V=orth(V);
    V(:,1)=ones(Nf,1);
    V = GramSchmidt(V);
    Ei=[1, 0.2.*ones(1,Nf-1)];
    Ei = sort(Ei,'descend');
    Di = diag(Ei);
    A = V*Di*V';
    A(A<0)=0;
    for ii = 1:Nf
        A(ii,:) =  A(ii,:)./sum(A(ii,:));
    end
    
    
    
    
    init_x(1,:) = X1{start_t};
    init_x(2,:) = X2{start_t};
    init_x(3,:) = X3{start_t};
    
    for ii = 1:Nf
        for kk = 1:Nm
            for jj = 1:Np
                PF.X{ii,kk}(jj,:) = mvnrnd(init_x(ii,:),[2,0,0;...
                    0,2,0;...
                    0,0,0.1]);
            end
        end
    end
    
    t = start_t;
    end_t = start_t+100;
    map_obs = [];
    while t<end_t
        
        
        Xt{1,2} =  X1{t}-X2{t};
        Xt{2,1} =  X2{t}-X1{t};
        Xt{1,3} =  X1{t}-X3{t};
        Xt{3,1} =  X3{t}-X1{t};
        Xt{2,3} =  X2{t}-X3{t};
        Xt{3,2} =  X3{t}-X2{t};
        Xt{1,1} =  X1{t}-X1{t};
        Xt{2,2} =  X1{t}-X1{t};
        Xt{3,3} =  X1{t}-X1{t};
        
        scandata1 = scans1{t};
        
        V(1) = omometry_data1{t-1}(1);
        
        omega(1) = omometry_data1{t-1}(2);
        n1 = length(scandata1.Ranges);
        currentScan{1}  = lidarScan(scans1{t});
        currentScan{1} = removeInvalidData(currentScan{1},'RangeLimits',[0.1 12]);
        
        scandata2 = scans2{t};
        V(2) = omometry_data2{t-1}(1);
        omega(2) = omometry_data2{t-1}(2);
        n2 = length(scandata2.Ranges);
        currentScan{2}  = lidarScan(scans2{t});
        currentScan{2} = removeInvalidData(currentScan{2},'RangeLimits',[0.1 12]);
        
        scandata3 = scans3{t};
        V(3) = omometry_data3{t-1}(1);
        omega(3) = omometry_data3{t-1}(2);
        n3 = length(scandata3.Ranges);
        currentScan{3}  = lidarScan(scans3{t});
        currentScan{3} = removeInvalidData(currentScan{3},'RangeLimits',[0.1 12]);
        
        %% model new step
        for ii = 1:Nf
            for kk = 1:Nm
                for jj = 1:Np
                    PF.X{ii,kk}(jj,:) = PF.X{ii,kk}(jj,:) + [(V(ii)+normrnd(0,sigma_x_step))*cos( PF.X{ii,kk}(jj,3)),(V(ii)+normrnd(0,sigma_y_step))*sin(PF.X{ii,kk}(jj,3)),omega(ii)].*dt{t}+mvnrnd(zeros(1,3),Cnoise);
                end
            end
        end
        %% Likelihood
        for ii = 1:Nf
            for kk = 1:Nm
                for qq = 1:Np
                        Relitive_scan = transformScan(currentScan{ii},PF.X{ii,kk}(qq,:));
                        estimated_obsticles = Relitive_scan.Cartesian;
                        PF.W{ii,kk}(qq) = Likelihood(estimated_obsticles,max_laser_samples,obsticle_vector,sigma_L);
                end
                if sum(isnan(PF.W{ii,kk}))>=1
                    PF.W{ii,kk} = ones(Np,1);
                end
                PF.W{ii,kk} = PF.W{ii,kk}./sum(PF.W{ii,kk});
            end
        end
        
        %% resmpling
        for ii = 1:Nf
            for kk = 1:Nm
                PF.X{ii,kk} = PF.X{ii,kk}(randsample(Np,Np,true,PF.W{ii,kk}),:);
                PF.W{ii,kk} = ones(Np,1)./Np;
            end
        end
        
        %% ------ Fusion using particles intersection ------ %%
        
        % ---- find alpha using trce ---- %
        a = trace_localization_V3(PF,1);
        w1_t = ParticlesIntersectionK_localization_V3( PF,1,a,0.02 ,Xt,1);
        w2_t = ParticlesIntersectionK_localization_V3( PF,1,a,0.02 ,Xt,2);
        w3_t = ParticlesIntersectionK_localization_V3( PF,1,a,0.02 ,Xt,3);
        PF.W{1,1} = w1_t./sum(w1_t);
        PF.W{2,1} = w2_t./sum(w2_t);
        PF.W{3,1} = w3_t./sum(w3_t);
        for ii = 1:Nf
            PF.X{ii,1} = PF.X{ii,1}(randsample(Np,Np,true,PF.W{ii,1}),:);
            PF.W{ii,1} = ones(Np,1)./Np;
        end
        
        
        
        
        
        % ---- find alpha using det ---- %
        a = det_localization_V3(PF,2);
        w1_t = ParticlesIntersectionK_localization_V3( PF,2,a,0.02 ,Xt,1);
        w2_t = ParticlesIntersectionK_localization_V3( PF,2,a,0.02 ,Xt,2);
        w3_t = ParticlesIntersectionK_localization_V3( PF,2,a,0.02 ,Xt,3);
        PF.W{1,2} = w1_t./sum(w1_t);
        PF.W{2,2} = w2_t./sum(w2_t);
        PF.W{3,2} = w3_t./sum(w3_t);
        for ii = 1:Nf
            PF.X{ii,2} = PF.X{ii,2}(randsample(Np,Np,true,PF.W{ii,2}),:);
            PF.W{ii,2} = ones(Np,1)./Np;
        end
        
        % ---- find alpha using extended det ---- %
        a = extended_det_localization_V3(PF,3);
        w1_t = ParticlesIntersectionK_localization_V3( PF,3,a,0.02 ,Xt,1);
        w2_t = ParticlesIntersectionK_localization_V3( PF,3,a,0.02 ,Xt,2);
        w3_t = ParticlesIntersectionK_localization_V3( PF,3,a,0.02 ,Xt,3);
        PF.W{1,3} = w1_t./sum(w1_t);
        PF.W{2,3} = w2_t./sum(w2_t);
        PF.W{3,3} = w3_t./sum(w3_t);
        for ii = 1:Nf
            PF.X{ii,3} = PF.X{ii,3}(randsample(Np,Np,true,PF.W{ii,3}),:);
            PF.W{ii,3} = ones(Np,1)./Np;
        end
        
        % ---- find alpha using maximun Chernoff information ---- %
        Chernoff_optim = @(x) -Chernoff_inf_localization_V3(PF,4,x,0.02,Xt,3,1000);
        x0 = [0.4; 0.3;0.3];
        A0 = eye(3);
        b0 = [1;1;1];
        Aeq = [1,1,1];
        beq = 1;
        lb = [0;0;0];
        ub = [1;1;1];
        [a,fval] = fmincon(Chernoff_optim,x0,A0,b0,Aeq,beq,lb,ub,[],options);
        w1_t = ParticlesIntersectionK_localization_V3( PF,4,a,0.02 ,Xt,1);
        w2_t = ParticlesIntersectionK_localization_V3( PF,4,a,0.02 ,Xt,2);
        w3_t = ParticlesIntersectionK_localization_V3( PF,4,a,0.02 ,Xt,3);
        PF.W{1,4} = w1_t./sum(w1_t);
        PF.W{2,4} = w2_t./sum(w2_t);
        PF.W{3,4} = w3_t./sum(w3_t);
        for ii = 1:Nf
            PF.X{ii,4} = PF.X{ii,4}(randsample(Np,Np,true,PF.W{ii,4}),:);
            PF.W{ii,4} = ones(Np,1)./Np;
        end
        
        
        
        %% ------ Fusion using Linear particles intersection ------ %%
        
        % ---- find alpha using trce ---- %
        a = trace_localization_V3(PF,5);
        x1_t = LinearParticlesIntersectionK_localization_V3( PF,5,a,0.02 ,Xt,1);
        x2_t = LinearParticlesIntersectionK_localization_V3( PF,5,a,0.02 ,Xt,2);
        x3_t = LinearParticlesIntersectionK_localization_V3( PF,5,a,0.02 ,Xt,3);
        PF.X{1,5} = x1_t;
        PF.X{2,5} = x2_t;
        PF.X{3,5} = x3_t;
        for ii = 1:Nf
            PF.W{ii,5} = ones(Np,1)./Np;
        end
        
        % ---- find alpha using det ---- %
        a = det_localization_V3(PF,6);
        x1_t = LinearParticlesIntersectionK_localization_V3( PF,6,a,0.02 ,Xt,1);
        x2_t = LinearParticlesIntersectionK_localization_V3( PF,6,a,0.02 ,Xt,2);
        x3_t = LinearParticlesIntersectionK_localization_V3( PF,6,a,0.02 ,Xt,3);
        PF.X{1,6} = x1_t;
        PF.X{2,6} = x2_t;
        PF.X{3,6} = x3_t;
        for ii = 1:Nf
            PF.W{ii,6} = ones(Np,1)./Np;
        end
        
        % ---- find alpha using extended det ---- %
        a = extended_det_localization_V3(PF,7);
        x1_t = LinearParticlesIntersectionK_localization_V3( PF,7,a,0.02 ,Xt,1);
        x2_t = LinearParticlesIntersectionK_localization_V3( PF,7,a,0.02 ,Xt,2);
        x3_t = LinearParticlesIntersectionK_localization_V3( PF,7,a,0.02 ,Xt,3);
        PF.X{1,7} = x1_t;
        PF.X{2,7} = x2_t;
        PF.X{3,7} = x3_t;
        for ii = 1:Nf
            PF.W{ii,7} = ones(Np,1)./Np;
        end
        
        % ---- find alpha using maximun Chernoff information ---- %
        Chernoff_optim = @(x) -Chernoff_inf_localization_V3(PF,8,x,0.02,Xt,3,1000);
        x0 = [0.4; 0.3;0.3];
        A0 = eye(3);
        b0 = [1;1;1];
        Aeq = [1,1,1];
        beq = 1;
        lb = [0;0;0];
        ub = [1;1;1];
        [a,fval] = fmincon(Chernoff_optim,x0,A0,b0,Aeq,beq,lb,ub,[],options);
        x1_t = LinearParticlesIntersectionK_localization_V3( PF,8,a,0.02 ,Xt,1);
        x2_t = LinearParticlesIntersectionK_localization_V3( PF,8,a,0.02 ,Xt,2);
        x3_t = LinearParticlesIntersectionK_localization_V3( PF,8,a,0.02 ,Xt,3);
        PF.X{1,8} = x1_t;
        PF.X{2,8} = x2_t;
        PF.X{3,8} = x3_t;
        for ii = 1:Nf
            PF.W{ii,8} = ones(Np,1)./Np;
        end
        
        
        
        %% ERROR
        XrTemp = [X1{t}(1:2);X2{t}(1:2);X3{t}(1:2)];
        for ii = 1:Nf
            for jj = 1:Nm
                PF_Error(t,Nmc,jj,ii) = norm(mean(PF.X{ii,jj}(:,1:2))-XrTemp(ii,:));
            end
        end
        
        
        
        %% PLOT visualization
        
        
        %
        %         figure(1)
        %         hold on;
        %         show(map);
        %         for ii = 1:Nf
        %             for jj = 1:Nm
        %                 plot(mean(PF.X{ii,jj}(:,1)),mean(PF.X{ii,jj}(:,2)),'+','MarkerSize',10);
        %             end
        %         end
        %         plot(X1{t}(1),X1{t}(2),'y+','MarkerSize',20);
        %         plot(X2{t}(1),X2{t}(2),'y+','MarkerSize',20);
        %         plot(X3{t}(1),X3{t}(2),'y+','MarkerSize',20);
        %         xlabel('X[m]');
        %         ylabel('Y[m]');
        %         xlim([-10.5 10.5]);
        %         ylim([-10.5 10.5]);
        %         hold off;
        %
        %
        %
        %         for ii = 1:Nf
        %             figure(ii+1)
        %             hold on;
        %             for jj = 1:Nm
        %                 plot(start_t:t,PF_Error(start_t:t,Nmc,jj,ii));
        %             end
        %             legend('PI trace', 'PI det', 'PI extended det', 'PI Chernoff',...
        %                 'LPI trace', 'LPI det', 'LPI extended det', 'LPI Chernoff','No fusion');
        %             hold off;
        %         end
        %
        %
        %
        %
        %
        %
        %         pause(0.0001);
        
        disp(['iter: ',num2str(Nmc),', step: ',num2str(t-start_t)]);
        t = t+1;
    end
    close all;
end
save(FileName,'PF_Error');

for ii = 1:Nf
    for jj = 1:Nm
        for kk = 1:100
            PF_plot(:,kk,jj,ii) = PF_Error(find(PF_Error(:,kk,jj,ii)~=0),kk,jj,ii);
        end
    end
end

for ii = 1:Nf
    figure(ii);
    hold on;
    for jj = 1:Nm
        plot(1:100,mean(PF_plot(:,:,jj,ii),2));
    end
    legend('PI trace', 'PI det', 'PI extended det', 'PI Chernoff',...
        'LPI trace', 'LPI det', 'LPI extended det', 'LPI Chernoff','No fusion')
    hold off;
end
figure(ii+1);
for ii = 1:Nf
    for jj = 1:Nm
        S(1) = std(mean(PF_plot(:,:,jj,ii),2));
        M(1) = mean(mean(PF_plot(:,:,jj,ii),2));
        bar((9*ii-(9-jj)),M,'FaceColor','none');
        hold on;
        h(jj) =  errorbar((9*ii-(9-jj)),M,S);
    end
    
end
legend([h],{'PI trace', 'PI det', 'PI extended det', 'PI Chernoff',...
    'LPI trace', 'LPI det', 'LPI extended det', 'LPI Chernoff','No fusion'});


markers = ['o','^','s','p','o','^','s','p','.'];
facecolor = ['r','r','r','r','b','b','b','b','none'];
for ii = 1:Nf
    figure(ii+4);
    for jj = 1:Nm
        S(1) = std(mean(PF_plot(:,:,jj,ii),2));
        M(1) = mean(mean(PF_plot(:,:,jj,ii),2));
        bar(jj,M,'FaceColor',facecolor(jj),'EdgeColor','black');
        hold on;
        h(jj) =  errorbar(jj,M,S,'Color','black','Marker',markers(jj),'MarkerSize',8,'MarkerFaceColor','white');
    end
    ylabel('MSE');
    set(gca,'XTick',[2.5 6.5 9],'XTickLabel',{'PI','LPI','No fusion'});
    legend([h],{'trace', 'det', 'extended det', 'Chernoff'});
end



