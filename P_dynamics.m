clear all;
diary log_cmdwindow
tStart=tic;

%% Generation of white Gaussian noise
sigma_s2=1e-4%1e2;%1e-2;%1e-6;% %system noise
sigma_o2=1e0%1e2;%1e-4;%1e-2;% %observation noise

%% Parameters =======================================================
r = 6;%10;% %内部自由度，潜在変数
n = 1000;
pmin = 1;
pinc = 1;
pmax = 20;
ps   = pmin:pinc:pmax;
% ps = 7;%[2 11];
% p=5;
% m=1;%
num_ave = 1;
num_ave2 = 1;
CNT = 0; % Counter
CNT2 = 0;
CNT3 = 0;

%% Preparation of output directories ================================
srcdir = ('src');
workdir   = ('../work');%('work_200803/Arand_Crand_2_r10_s001_o1e-1_ave10');
workdir2  = [workdir,'/fig_loop'];
% videodir  = [workdir,'/video'];
% sensordir = [workdir,'/sensor_location'];
mkdir(workdir);
mkdir(workdir2);
% mkdir(videodir);
% mkdir(sensordir);

% %% Generation dynamical system ======================================
% % for p = ps
%     B_sys = 0;
%     while B_sys == 0
%         sysorg=drss(r,n); % generate random discrete test model
%                             % (r-th order model with one input and n outputs)
%         B_sys = isstable(sysorg); % Boolean value of stability
%         CNT2=CNT2+1;
%     end
%     [Aorg,Borg,Corg,Dorg]=ssdata(sysorg); % extract the matrix data A,B,C,D from the state-model sys
%     % [Atmp,Borg,Corg,Dorg]=ssdata(sysorg); 
%     % text = ['work_200720/A1_C1_s001_o001/Aorg.mat'];%[workdir,'/Arog.mat'];
%     % load(text);
%     % text = ['work_200720/A1_C1_s001_o001/Corg.mat'];%[workdir,'/Arog.mat'];
%     % load(text);
%     % Aeigmax = 1e1
%     % while Aeigmax > 1
%     %     sysorg=rss(r,n);
%     %     [Aorg,Borg,Corg,Dorg]=ssdata(sysorg);
%     %     Aeigmax = max(eig(Aorg))
%     % end
% % end

%% Initial conditions ===============================================
% H=[];
% H_2=[];
% X=zeros(r,m);
% X(:,1)=randn(r,1);%初期時刻は与えるよね？
% Y=zeros(p,m); %使わない？
% B=eye(r,1); %カルマンフィルタではこうなる？
% D=[]; %？

% noise=wgn(n,1,0) %noise=wgn(m,n,power): m行n列の行列・dB W単位のパワー指定？

%% Average loop =====================================================
rng(1,'twister');
for w=1:1:num_ave2
    % if mod(w,num_ave2/10) == 0
    %     disp(['Average loop: ',num2str(w),'/',num2str(num_ave2)]);
    % end
    text = [ 'Average loop: ', num2str(w),'/', num2str(num_ave2) ]; %text = [ 'Average loop: ', num2str(w),' ...' ];
    disp(text);
    %% Generation dynamical system ==================================
    % % random system ///
    % B_sys = 0;
    % cond_AAI = 1e10;
    % CNT2=0;
    % while (B_sys == 0) || (cond_AAI > 1e3)
    %     CNT2 = CNT2+1;
    %     sysorg = drss(r,n); % generate random discrete test model
    %                         % (r-th order model with one input and n outputs)
    %     B_sys = isstable(sysorg); % Boolean value of stability
    %     [Aorg,Borg,Corg,Dorg] = ssdata(sysorg);
    %     AAI = kron(Aorg',Aorg')-eye(r^2);
    %     cond_AAI = cond(AAI);
    %     text = [ 'CNT2=', num2str(CNT2), ', B_sys=', num2str(B_sys),...
    %              ', cond(AAI)=', num2str(cond_AAI) ]; 
    %     disp(text);
    % end

    % linear system of three dynamic modes ///
    % specify characteristic frequencies and growth/decay rates
    % associated with continuous-time dynamics
    f = [1.0 2.5 5.5];
%    g = [-0.0001 -0.0001 -0.003];%[0 0 -.3];%[0.001 0.001 0.03];%[0 0 0.3];%
    g = [-50.0 -50.0 -1500.0];
    g = [-100.0 -100.0 -3000.0];
    assert(length(f)==length(g))
    % construct low-rank continuous-time operator (rank=k)
    k = 2*length(f);  % (to perform DMD/TDMD with correct rank, set r=k)
    A1 = [];
    for ii = 1:length(f)
        A2 = [[g(ii) 2*pi*f(ii); -2*pi*f(ii) g(ii)]];
        A1 = [A1 A2];
    end
    Alowrank = [];
    for ii = 1:length(f)
        Alowrank = blkdiag(Alowrank,A1(:,(ii-1)*2+1:2*ii));
    end
    % size(Alowrank)
    dt = 0.01; % time step size
    ii = 1;
    Aorg = expm(dt*ii*Alowrank);
    AAI = kron(Aorg',Aorg')-eye(r^2);
    sysorg = drss(r,n);
    [~,Borg,Corg,Dorg] = ssdata(sysorg);

    % ///
    sysmemo(w,1) = norm(Aorg(:));
    sysmemo(w,2) = norm(Borg(:));
    sysmemo(w,3) = norm(Corg(:));
    sysmemo(w,4) = norm(Dorg(:));

    eigA = eig(Aorg)
    absA = abs(eigA)
    
    %% Sensor selection =============================================
    CNT=0;
    for p = ps
        CNT = CNT+1;
        text = [ num2str(p),' sensor selection started --->' ]; 
        disp(text);

        % %% Average loop =============================================
        % for w=1:1:num_ave

        %     %% Generation dynamical system ======================================
        %     B_sys = 0;
        %     while B_sys == 0
        %         sysorg = drss(r,n); % generate random discrete test model
        %         B_sys = isstable(sysorg); %Boolean value of stability
        %         % Aeigmax = max(eig(Aorg));
        %         % text = [ 'Aeigmax:', num2str(Aeigmax), 'sys:' , num2str(B_sys) ];
        %         % disp(text);
        %     end
        %     [Aorg,Borg,Corg,Dorg]=ssdata(sysorg);
        %     % [~,~,Corg,~]=ssdata(sysorg); %remake only Corg
        %     % sysorg_check = ss(Aorg,Borg,Corg,Dorg);
        %     % B_sys = isstable(sysorg_check)

            %% Dynamics w/ noise-----------------------------------------
            % Minimization of steady state error covariance matrix (Kalman filter)
            [time_KF(CNT,w+1), H_KF, sensors_KF] = F_sensor_KF(Aorg,Borg,Corg,Dorg,p,sigma_s2,sigma_o2);
            detWo_KF   (CNT,w+1) = F_calc_detWo (p,H_KF,Aorg,Borg,Corg,Dorg);
            detWo_KF_2 (CNT,w+1) = F_calc_detWo2(p,H_KF,Aorg,Borg,Corg,Dorg);
            detWo_KF_3 (CNT,w+1) = F_calc_detWo3(p,H_KF,Aorg,Borg,Corg,Dorg);
            detP_KF (CNT,w+1) = F_calc_detP  (p,H_KF,Aorg,Borg,Corg,Dorg,sigma_s2,sigma_o2);
            detCC_KF (CNT,w+1) = F_calc_det (p,H_KF,Corg);

            %% Dynamics w/o noise -----------------------------------------
            % Maximization of obserbility Gramian
            [time_Gram(CNT,w+1), H_Gram, sensors_Gram] = F_sensor_Gram(Aorg,Borg,Corg,Dorg,p,sigma_s2,sigma_o2);
            detWo_Gram   (CNT,w+1) = F_calc_detWo (p,H_Gram,Aorg,Borg,Corg,Dorg);
            detWo_Gram_2 (CNT,w+1) = F_calc_detWo2(p,H_Gram,Aorg,Borg,Corg,Dorg);
            detWo_Gram_3 (CNT,w+1) = F_calc_detWo3(p,H_Gram,Aorg,Borg,Corg,Dorg);
            detP_Gram (CNT,w+1) = F_calc_detP  (p,H_Gram,Aorg,Borg,Corg,Dorg,sigma_s2,sigma_o2);
            detCC_Gram (CNT,w+1) = F_calc_det (p,H_Gram,Corg);
            
            %% S2S w/ noise-----------------------------------------
            % Maximization of determinant of (CRC)
            [time_S2SwR(CNT,w+1), H_S2SwR, sensors_S2SwR] = F_sensor_DGwR(Corg,p,sigma_o2);
            detWo_S2SwR   (CNT,w+1) = F_calc_detWo (p,H_S2SwR,Aorg,Borg,Corg,Dorg);
            detWo_S2SwR_2 (CNT,w+1) = F_calc_detWo2(p,H_S2SwR,Aorg,Borg,Corg,Dorg);
            detWo_S2SwR_3 (CNT,w+1) = F_calc_detWo3(p,H_S2SwR,Aorg,Borg,Corg,Dorg);
            detP_S2SwR (CNT,w+1) = F_calc_detP  (p,H_S2SwR,Aorg,Borg,Corg,Dorg,sigma_s2,sigma_o2);
            detCC_S2SwR (CNT,w+1) = F_calc_det (p,H_S2SwR,Corg);
            
            %% S2S w/o noise-----------------------------------------
            % Maximization of determinant of Fisher information matrix
            [time_S2S(CNT,w+1), H_S2S, sensors_S2S] = F_sensor_QD(Corg,p);
            detWo_S2S   (CNT,w+1) = F_calc_detWo (p,H_S2S,Aorg,Borg,Corg,Dorg);
            detWo_S2S_2 (CNT,w+1) = F_calc_detWo2(p,H_S2S,Aorg,Borg,Corg,Dorg);
            detWo_S2S_3 (CNT,w+1) = F_calc_detWo3(p,H_S2S,Aorg,Borg,Corg,Dorg);
            detP_S2S (CNT,w+1) = F_calc_detP  (p,H_S2S,Aorg,Borg,Corg,Dorg,sigma_s2,sigma_o2);
            detCC_S2S (CNT,w+1) = F_calc_det (p,H_S2S,Corg);

        % end
        
        % %% Averaging ================================================
        % [ time_KF, detWo_KF, detP_KF ]...
        % = F_data_ave1( CNT, num_ave, time_KF, detWo_KF, detP_KF );
        % [ time_std_KF, detWo_std_KF, detP_std_KF ]...
        % = F_data_std1( CNT, num_ave, time_KF, detWo_KF, detP_KF );
        % [ time_Gram, detWo_Gram, detP_Gram ]...
        % = F_data_ave1( CNT, num_ave, time_Gram, detWo_Gram, detP_Gram );
        % [ time_std_Gram, detWo_std_Gram, detP_std_Gram ]...
        % = F_data_std1( CNT, num_ave, time_Gram, detWo_Gram, detP_Gram );
        % [ time_S2SwR, detWo_S2SwR, detP_S2SwR ]...
        % = F_data_ave1( CNT, num_ave, time_S2SwR, detWo_S2SwR, detP_S2SwR );
        % [ time_std_S2SwR, detWo_std_S2SwR, detP_std_S2SwR ]...
        % = F_data_std1( CNT, num_ave, time_S2SwR, detWo_S2SwR, detP_S2SwR );
        % [ time_S2S, detWo_S2S, detP_S2S ]...
        % = F_data_ave1( CNT, num_ave, time_S2S, detWo_S2S, detP_S2S );
        % [ time_std_S2S, detWo_std_S2S, detP_std_S2S ]...
        % = F_data_std1( CNT, num_ave, time_S2S, detWo_S2S, detP_S2S );

        text = [ '---> ', num2str(p), ' sensor selection finished!' ];
        % disp(text);
    end

    %% Plot =====================================================
    cd(workdir2)
    newcolors = {'red','blue','green','black'};
    colororder(newcolors)
    fig1 = semilogy(ps,detP_KF(1:CNT,w+1),'-o','DisplayName','KF-based-greedy','LineWidth',1,'MarkerSize',5);
    hold on
    fig1 = semilogy(ps,detP_Gram(1:CNT,w+1),'-s','DisplayName','Gramian-based-greedy','LineWidth',1,'MarkerSize',5);
    hold on
    % fig1 = semilogy(ps,detP_S2SwR(1:CNT,w+1),'-*','DisplayName','S2SwR','LineWidth',1,'MarkerSize',5);
    % hold on
    fig1 = semilogy(ps,detP_S2S(1:CNT,w+1),'-+','DisplayName','S2S-based-greedy','LineWidth',1,'MarkerSize',5);
    hold off
    legend({}, 'FontSize', 14)
    ax = gca;
    ax.FontSize = 14;
    xlabel('Number of sensors', 'FontSize', 14)
    ylabel('det P', 'FontSize', 14)
    filename = [ 'loop', num2str(w), '_detP.png' ];
    saveas(fig1, filename);
    
    colororder(newcolors)
    fig2 = semilogy(ps,detWo_KF(1:CNT,w+1),'-o','DisplayName','KF-based-greedy','LineWidth',1,'MarkerSize',5);
    hold on
    fig2 = semilogy(ps,detWo_Gram(1:CNT,w+1),'-s','DisplayName','Gramian-based-greedy','LineWidth',1,'MarkerSize',5);
    hold on
    % fig1 = semilogy(ps,detWo_S2SwR(1:CNT,w+1),'-*','DisplayName','S2SwR','LineWidth',1,'MarkerSize',5);
    % hold on
    fig2 = semilogy(ps,detWo_S2S(1:CNT,w+1),'-+','DisplayName','S2S-based-greedy','LineWidth',1,'MarkerSize',5);
    hold off
    legend({}, 'FontSize', 14, 'Location', 'southeast')
    ax = gca;
    ax.FontSize = 14;
    xlabel('Number of sensors', 'FontSize', 14)
    ylabel('det W_o', 'FontSize', 14)
    filename = [ 'loop', num2str(w), '_detWo.png' ];
    saveas(fig2, filename);

    colororder(newcolors)
    fig3 = semilogy(ps,detCC_KF(1:CNT,w+1),'-o','DisplayName','KF-based-greedy','LineWidth',1,'MarkerSize',5);
    hold on
    fig3 = semilogy(ps,detCC_Gram(1:CNT,w+1),'-s','DisplayName','Gramian-based-greedy','LineWidth',1,'MarkerSize',5);
    hold on
    % fig1 = semilogy(ps,detCC_S2SwR(1:CNT,w+1),'-*','DisplayName','S2SwR','LineWidth',1,'MarkerSize',5);
    % hold on
    fig3 = semilogy(ps,detCC_S2S(1:CNT,w+1),'-+','DisplayName','S2S-based-greedy','LineWidth',1,'MarkerSize',5);
    hold off
    legend({}, 'FontSize', 14, 'Location', 'southeast')
    ax = gca;
    ax.FontSize = 14;
    xlabel('Number of sensors', 'FontSize', 14)
    ylabel('det CC^T or det C^TC', 'FontSize', 14)
    filename = [ 'loop', num2str(w), '_detCC.png' ];
    saveas(fig3, filename);
    
    colororder(newcolors)
    fig4 = semilogy(ps,1./detP_KF(1:CNT,w+1),'-o','DisplayName','KF-based-greedy','LineWidth',1,'MarkerSize',5);
    hold on
    fig4 = semilogy(ps,1./detP_Gram(1:CNT,w+1),'-s','DisplayName','Gramian-based-greedy','LineWidth',1,'MarkerSize',5);
    hold on
    % fig1 = semilogy(ps,1./detP_S2SwR(1:CNT,w+1),'-*','DisplayName','S2SwR','LineWidth',1,'MarkerSize',5);
    % hold on
    fig4 = semilogy(ps,1./detP_S2S(1:CNT,w+1),'-+','DisplayName','S2S-based-greedy','LineWidth',1,'MarkerSize',5);
    hold off
    legend({}, 'FontSize', 14, 'Location', 'southeast')
    ax = gca;
    ax.FontSize = 14;
    xlabel('Number of sensors', 'FontSize', 14)
    ylabel('det P^{-1}', 'FontSize', 14)
    filename = [ 'loop', num2str(w), '_detP-1.png' ];
    saveas(fig4, filename);
    % cd ../../
    cd ../
    save('sysmemo.mat','sysmemo')
    cd ../
    save('eig.mat','eigA','absA')
    cd ../
    cd(srcdir)

end

%% Averaging ================================================
[ time_KF, detWo_KF, detWo_KF_2, detWo_KF_3, detP_KF, detCC_KF ]...
= F_data_ave3( CNT, num_ave2, time_KF, detWo_KF, detWo_KF_2, detWo_KF_3, detP_KF, detCC_KF );
[ time_std_KF, detWo_std_KF, detWo_std_KF_2, detWo_std_KF_3, detP_std_KF, detCC_std_KF ]...
= F_data_std2( CNT, num_ave2, time_KF, detWo_KF, detWo_KF_2, detWo_KF_3, detP_KF, detCC_KF );
[ time_Gram, detWo_Gram, detWo_Gram_2, detWo_Gram_3, detP_Gram, detCC_Gram ]...
= F_data_ave3( CNT, num_ave2, time_Gram, detWo_Gram, detWo_Gram_2, detWo_Gram_3, detP_Gram, detCC_Gram );
[ time_std_Gram, detWo_std_Gram, detWo_std_Gram_2, detWo_std_Gram_3, detP_std_Gram, detCC_std_Gram ]...
= F_data_std2( CNT, num_ave2, time_Gram, detWo_Gram, detWo_Gram_2, detWo_Gram_3, detP_Gram, detCC_Gram );
[ time_S2SwR, detWo_S2SwR, detWo_S2SwR_2, detWo_S2SwR_3, detP_S2SwR, detCC_S2SwR ]...
= F_data_ave3( CNT, num_ave2, time_S2SwR, detWo_S2SwR, detWo_S2SwR_2, detWo_S2SwR_3, detP_S2SwR, detCC_S2SwR );
[ time_std_S2SwR, detWo_std_S2SwR, detWo_std_S2SwR_2, detWo_std_S2SwR_3, detP_std_S2SwR, detCC_std_S2SwR ]...
= F_data_std2( CNT, num_ave2, time_S2SwR, detWo_S2SwR, detWo_S2SwR_2, detWo_S2SwR_3, detP_S2SwR, detCC_S2SwR );
[ time_S2S, detWo_S2S, detWo_S2S_2, detWo_S2S_3, detP_S2S, detCC_S2S ]...
= F_data_ave3( CNT, num_ave2, time_S2S, detWo_S2S, detWo_S2S_2, detWo_S2S_3, detP_S2S, detCC_S2S );
[ time_std_S2S, detWo_std_S2S, detWo_std_S2S_2, detWo_std_S2S_3, detP_std_S2S, detCC_std_S2S ]...
= F_data_std2( CNT, num_ave2, time_S2S, detWo_S2S, detWo_S2S_2, detWo_S2S_3, detP_S2S, detCC_S2S );

% Data organization =================================================
% Arrange
[time_all]  = F_data_arrange1( ps, CNT, time_KF, time_Gram, time_S2SwR, time_S2S );
[detWo_all] = F_data_arrange1( ps, CNT, detWo_KF, detWo_Gram, detWo_S2SwR, detWo_S2S );
[detWo2_all]= F_data_arrange1( ps, CNT, detWo_KF_2, detWo_Gram_2, detWo_S2SwR_2, detWo_S2S_2 );
[detWo3_all]= F_data_arrange1( ps, CNT, detWo_KF_3, detWo_Gram_3, detWo_S2SwR_3, detWo_S2S_3 );
[detP_all]  = F_data_arrange1( ps, CNT, detP_KF, detP_Gram, detP_S2SwR, detP_S2S );
[detCC_all] = F_data_arrange1( ps, CNT, detCC_KF, detCC_Gram, detCC_S2SwR, detCC_S2S );
[time_std]  = F_data_arrange1( ps, CNT, time_std_KF, time_std_Gram, time_std_S2SwR, time_std_S2S );
[detWo_std] = F_data_arrange1( ps, CNT, detWo_std_KF, detWo_std_Gram, detWo_std_S2SwR, detWo_std_S2S );
[detWo2_std]= F_data_arrange1( ps, CNT, detWo_std_KF_2, detWo_std_Gram_2, detWo_std_S2SwR_2, detWo_std_S2S_2 );
[detWo3_std]= F_data_arrange1( ps, CNT, detWo_std_KF_3, detWo_std_Gram_3, detWo_std_S2SwR_3, detWo_std_S2S_3 );
[detP_std]  = F_data_arrange1( ps, CNT, detP_std_KF, detP_std_Gram, detP_std_S2SwR, detP_std_S2S );
[detCC_std] = F_data_arrange1( ps, CNT, detCC_std_KF, detCC_std_Gram, detCC_std_S2SwR, detCC_std_S2S );

% Normalize


%% Save =============================================================
cd(workdir)
% Average
save('time.mat','time_all');
save('detWo.mat','detWo_all');
save('detWo2.mat','detWo2_all');
save('detWo3.mat','detWo3_all');
save('detP.mat','detP_all');
save('detCC.mat','detCC_all');
% std
save('time_std.mat','time_std');
save('detWo_std.mat','detWo_std');
save('detWo2_std.mat','detWo2_std');
save('detWo3_std.mat','detWo3_std');
save('detP_std.mat','detP_std');
save('detCC_std.mat','detCC_std');
% each detP
save('detP_KF.mat','detP_KF');
save('detP_Gram.mat','detP_Gram');
save('detP_S2SwR.mat','detP_S2SwR');
save('detP_S2S.mat','detP_S2S');
% each detWo
save('detWo_KF.mat','detWo_KF');
save('detWo_Gram.mat','detWo_Gram');
save('detWo_S2SwR.mat','detWo_S2SwR');
save('detWo_S2S.mat','detWo_S2S');
% each detWo2
save('detWo_KF_2.mat','detWo_KF_2');
save('detWo_Gram_2.mat','detWo_Gram_2');
save('detWo_S2SwR_2.mat','detWo_S2SwR_2');
save('detWo_S2S_2.mat','detWo_S2S_2');
% each detWo3
save('detWo_KF_3.mat','detWo_KF_3');
save('detWo_Gram_3.mat','detWo_Gram_3');
save('detWo_S2SwR_3.mat','detWo_S2SwR_3');
save('detWo_S2S_3.mat','detWo_S2S_3');
% each detCC
save('detCC_KF.mat','detCC_KF');
save('detCC_Gram.mat','detCC_Gram');
save('detCC_S2SwR.mat','detCC_S2SwR');
save('detCC_S2S.mat','detCC_S2S');
% each time
save('time_KF','time_KF');
save('time_Gram','time_Gram');
save('time_S2SwR','time_S2SwR');
save('time_S2S','time_S2S');

%% Plot =============================================================
newcolors = {'red','blue','green','black'};
colororder(newcolors)
fig1 = semilogy(ps,detP_all(1:CNT,2),'-o','DisplayName','KF-based-greedy','LineWidth',1,'MarkerSize',5);
hold on
fig1 = semilogy(ps,detP_all(1:CNT,3),'-s','DisplayName','Gramian-based-greedy','LineWidth',1,'MarkerSize',5);
hold on
% fig1 = semilogy(ps,detP_all(1:CNT,4),'-*','DisplayName','S2SwR','LineWidth',1,'MarkerSize',5);
% hold on
fig1 = semilogy(ps,detP_all(1:CNT,5),'-+','DisplayName','S2S-based-greedy','LineWidth',1,'MarkerSize',5);
hold off
legend({}, 'FontSize', 14)
ax = gca;
ax.FontSize = 14;
xlabel('Number of sensors', 'FontSize', 14)
ylabel('det P', 'FontSize', 14)
saveas(fig1,'fig_detP.png');
    
colororder(newcolors)
fig7 = semilogy(ps,1./detP_all(1:CNT,2),'-o','DisplayName','KF-based-greedy','LineWidth',1,'MarkerSize',5);
hold on
fig7 = semilogy(ps,1./detP_all(1:CNT,3),'-s','DisplayName','Gramian-based-greedy','LineWidth',1,'MarkerSize',5);
hold on
fig7 = semilogy(ps,1./detP_all(1:CNT,5),'-+','DisplayName','S2S-based-greedy','LineWidth',1,'MarkerSize',5);
hold off
legend({}, 'FontSize', 14, 'Location', 'southeast')
ax = gca;
ax.FontSize = 14;
xlabel('Number of sensors', 'FontSize', 14)
ylabel('det P^{-1}', 'FontSize', 14)
filename = [ 'fig_detP-1.png' ];
saveas(fig7, filename);

colororder(newcolors)
fig2 = semilogy(ps,detWo_all(1:CNT,2),'-o','DisplayName','KF-based-greedy','LineWidth',1,'MarkerSize',5);
hold on
fig2 = semilogy(ps,detWo_all(1:CNT,3),'-s','DisplayName','Gramian-based-greedy','LineWidth',1,'MarkerSize',5);
hold on
% fig2 = semilogy(ps,detWo_all(1:CNT,4),'-*','DisplayName','S2SwR','LineWidth',1,'MarkerSize',5);
% hold on
fig2 = semilogy(ps,detWo_all(1:CNT,5),'-+','DisplayName','S2S-based-greedy','LineWidth',1,'MarkerSize',5);
hold off
legend({}, 'FontSize', 14, 'Location', 'southeast')
ax = gca;
ax.FontSize = 14;
xlabel('Number of sensors', 'FontSize', 14)
ylabel('det W_o', 'FontSize', 14)
saveas(fig2,'fig_detWo.png');

colororder(newcolors)
fig3 = semilogy(ps,detWo2_all(1:CNT,2),'-o','DisplayName','KF-based-greedy','LineWidth',1,'MarkerSize',5);
hold on
fig3 = semilogy(ps,detWo2_all(1:CNT,3),'-s','DisplayName','Gramian-based-greedy','LineWidth',1,'MarkerSize',5);
hold on
% fig2 = semilogy(ps,detWo2_all(1:CNT,4),'-*','DisplayName','S2SwR','LineWidth',1,'MarkerSize',5);
% hold on
fig3 = semilogy(ps,detWo2_all(1:CNT,5),'-+','DisplayName','S2S-based-greedy','LineWidth',1,'MarkerSize',5);
hold off
legend({}, 'FontSize', 14, 'Location', 'southeast')
ax = gca;
ax.FontSize = 14;
xlabel('Number of sensors', 'FontSize', 14)
ylabel('det W_o', 'FontSize', 14)
saveas(fig3,'fig_detWo2.png');

colororder(newcolors)
fig4 = semilogy(ps,detWo3_all(1:CNT,2),'-o','DisplayName','KF-based-greedy','LineWidth',1,'MarkerSize',5);
hold on
fig4 = semilogy(ps,detWo3_all(1:CNT,3),'-s','DisplayName','Gramian-based-greedy','LineWidth',1,'MarkerSize',5);
hold on
% fig2 = semilogy(ps,detWo3_all(1:CNT,4),'-*','DisplayName','S2SwR','LineWidth',1,'MarkerSize',5);
% hold on
fig4 = semilogy(ps,detWo3_all(1:CNT,5),'-+','DisplayName','S2S-based-greedy','LineWidth',1,'MarkerSize',5);
hold off
legend({}, 'FontSize', 14, 'Location', 'southeast')
ax = gca;
ax.FontSize = 14;
xlabel('Number of sensors', 'FontSize', 14)
ylabel('det W_o', 'FontSize', 14)
saveas(fig4,'fig_detWo3.png');

newcolors = {'red','blue','green','black'};
colororder(newcolors)
fig5 = semilogy(ps,detCC_all(1:CNT,2),'-o','DisplayName','KF-based-greedy','LineWidth',1,'MarkerSize',5);
hold on
fig5 = semilogy(ps,detCC_all(1:CNT,3),'-s','DisplayName','Gramian-based-greedy','LineWidth',1,'MarkerSize',5);
hold on
% fig1 = semilogy(ps,detP_all(1:CNT,4),'-*','DisplayName','S2SwR','LineWidth',1,'MarkerSize',5);
% hold on
fig5 = semilogy(ps,detCC_all(1:CNT,5),'-+','DisplayName','S2S-based-greedy','LineWidth',1,'MarkerSize',5);
hold off
legend({}, 'FontSize', 14, 'Location', 'southeast')
ax = gca;
ax.FontSize = 14;
xlabel('Number of sensors', 'FontSize', 14)
ylabel('det CC^T or det C^TC', 'FontSize', 14)
saveas(fig5,'fig_detCC.png');

colororder(newcolors)
fig6 = semilogy(ps,time_all(1:CNT,2),'-o','DisplayName','KF-based-greedy','LineWidth',1,'MarkerSize',5);
hold on
fig6 = semilogy(ps,time_all(1:CNT,3),'-s','DisplayName','Gramian-based-greedy','LineWidth',1,'MarkerSize',5);
hold on
% fig3 = semilogy(ps,time_all(1:CNT,4),'-*','DisplayName','S2SwR','LineWidth',1,'MarkerSize',5);
% hold on
fig6 = semilogy(ps,time_all(1:CNT,5),'-+','DisplayName','S2S-based-greedy','LineWidth',1,'MarkerSize',5);
hold off
legend({}, 'FontSize', 14, 'Location', 'southeast')
ax = gca;
ax.FontSize = 14;
xlabel('Number of sensors', 'FontSize', 14)
ylabel('Computational time', 'FontSize', 14)
saveas(fig6,'fig_time.png');

disp('Congratulations!');
diary off
tEnd=toc(tStart)
save('runtime.mat','tEnd')
cd ../src%../../%
%% ///////////////////////////////////////////////////////////////////
%% Main program end


