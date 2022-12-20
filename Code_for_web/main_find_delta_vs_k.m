clear; format compact;
% close all; clc;
fold = 'c:\Users\user\Dropbox\Havlin\DynamicsRecovering/single_node_recovery/code\';
addpath([fold,'functions']);

global A Anw x0 Factor_Xeff List lambda NameOfModel a b

tic

% model of dynamics
NameOfModel = 'MM'; % 'Eco'; % 'Simple'; % 'Glauber'; % 'Voter'; % 'SIS';% 'MAK';%  'PD';%

M = KindOfDynamics( NameOfModel );
range = [0,100];

ChooseExcNodes = 'rand'; % options - {'low degree', 'high degree', 'local', 'rand','1'}
conditions.type = 'BC';% 'IC'; %

conds = {'free high','ignite'};
release = 1;
all_or_eff = 'eff'; % what xss we need: 'all'/ 'eff' /'eff_F'/ 'average'

x_th = 1/2;
NVec = 1e4; %round(logspace(1,4,11)); %1e5*(1:10);
lambdaVec = 3; %linspace(2.2,3,5);

w = 0.6;
k0Vec = logspace(log10(1.2),log10(40),20);
lmax = 3; % only for tree
kRealVec = zeros(size(k0Vec));
k2 = kRealVec; % second moment

holding_value = 2;
free_value_vec = [0 10];
NExcVec = [0 1]; %unique( round(logspace(0,log10(NVec/25),21)) );
DeltaCMat = zeros(1,length(k0Vec));
DeltaCMatstd = DeltaCMat;
Critical_delta = 1;
reals = 10;
CrDelRealVec = zeros(1,reals);

for idxN = 1:length(NVec)
    N = NVec(idxN)
    for idxk = 1:length(k0Vec)
        k0 = k0Vec(idxk)
        % build the network
        for idxl = 1:length(lambdaVec)
            NetStruct = 'ER';
            parameters = k0;
            if strcmp(NetStruct,'SF')
                lambda = lambdaVec(idxl);
                parameters = [lambda k0];
            end
            
            
            for idx_cond = 1:length(conds)
                    cond = conds{idx_cond};
                    switch cond
                        case 'free high'
                            conditions.free_value = free_value_vec(2);
                            NumExcited = NExcVec(1);
                        case 'ignite'
                            conditions.free_value = free_value_vec(1);
                            NumExcited = NExcVec(2);
                    end
                
                for ir = 1:reals
                    if strcmp(NetStruct,'RT'); N=1+k0*((k0-1)^lmax - 1)/(k0-2); end
                    Anw = BuildNetwork(N, NetStruct, parameters,'gcc'); % Adjacency matrix, not weighted
                    n = size(Anw,1);
                    A = w*Anw;
                    Factor_Xeff = sum(A,1)/sum(sum(A));
                    
                    if strcmp(NetStruct,'SF')
                        kRealVec(idxl) = mean(sum(Anw)); k2(idxl) = mean(sum(Anw).^2);
                    else
                        kRealVec(idxk) = mean(sum(Anw));
                        k2(idxk) = mean(sum(Anw).^2);
                    end
                    x0 = Initial_or_fixed_Conditions(n,NumExcited,ChooseExcNodes,conditions,holding_value);
                    toplot = 0;
                    
                    % find critical omega
                    options = optimset('TolX',1e-3);
                    % check limits
                    xA = find_x_via_delta(1e-2,M,all_or_eff,toplot,conditions,1);
                    xB = find_x_via_delta(50,M,all_or_eff,toplot,conditions,1);
                    if xB<x_th
                        Critical_delta = inf;
                    elseif xA>x_th
                        Critical_delta = nan;                        
                    else
                        Critical_delta = ...
                            fzero(@(delta) find_x_via_delta(delta,M,all_or_eff,toplot,conditions,1) - x_th ,[1e-3,50] ,options);
                    end
                    % PlotExamples when ...
                    for Examples=[]
                        if  (idxNE==1) && (ir==1) && (idxl==1) && (idxN==length(NVec)) && (idxk==1)
                            AxXvec = linspace(eps,2*Critical_delta,50);
                            AxYvec = arrayfun(@(omega) find_x_via_omega(omega,x0,M,all_or_eff,toplot,conditions,release),AxXvec);
                            arrayfun(@(omega) find_x_via_omega(omega,x0,M,'eff',1,conditions,release),Critical_delta*[0.99,1.01]);
                            
                            figure; % Xeff vs omega
                            plot(AxXvec, AxYvec,'.'...
                                ,AxXvec,(Factor_Xeff*x0)*ones(size(AxXvec)),'.')
                            xlabel('\omega_c'); ylabel('x_{eff}'); legend('SteadyState','InitialCondition');
                            title(['k = ',num2str(kRealVec(idxk)),', \rho = ',num2str(num2str(NExcVec(idxNE)/n))]);
                            set(gca,'XTick',sort([xlim,Critical_delta]),'XGrid','on','box','on','fontsize',15);
                            
                            figure; % Xeff vs beta
                            plot(AxXvec*k2(idxk)/kRealVec(idxk),AxYvec,'.'...
                                ,AxXvec*k2(idxk)/kRealVec(idxk),(Factor_Xeff*x0)*ones(size(AxXvec)),'.')
                            xlabel('\beta_c'); ylabel('x_{eff}'); legend('SteadyState','InitialCondition');
                            title(['k = ',num2str(kRealVec(idxk)),', \rho = ',num2str(num2str(NExcVec(idxNE)/n))]);
                            set(gca,'fontsize',15); box on
                        end
                    end
                    
                    % X0eff = Factor_Xeff*x0;
                    % beta = omega*sum(sum(A^2))/sum(A(:));
                    CrDelRealVec(ir) = Critical_delta;
                end
                if strcmp(NetStruct,'SF')
                    DeltaCMat(idx_cond, idxl) = mean(CrDelRealVec);
                    DeltaCMatstd(idx_cond,idxl) = std(CrDelRealVec);
                else
                    DeltaCMat(idx_cond,idxk) = mean(CrDelRealVec);
                    DeltaCMatstd(idx_cond,idxk) = std(CrDelRealVec);
                end
            end
        end
    end
    toc
end

kappa = k2./kRealVec-1;

dc_theory = arrayfun(@(kappa) find_dc_by_wk_theory_intersection(w,kappa,range) ,  kappa);

%% 

figure; hold on

set(gca,'FontSize',20,'box','on','XTick',[10 20 30],'LineWidth',1.5,'layer','top','XScale','log','YScale','log')
xlabel('\boldmath$\kappa$','Interpreter','latex','FontSize',33)
ylabel('\boldmath$ \Delta$','Interpreter','latex','FontSize',33)
axis square
xlim(kappa([1,end]))
ymin = 1e-1; ymax = 10^1.5;
ylim([ymin ymax])

% patches
clrs = [0 51 102; 255 192 0; 155 0 0]/255;

kappa_no_nan = kappa(isfinite(DeltaCMat(2,:)));
delte_no_nan = DeltaCMat(2,isfinite(DeltaCMat(2,:)));
idxc_inactive = find(isnan(DeltaCMat(1,:)),1);

patch([kappa,fliplr(kappa)],[ymin+0*kappa,ymax+0*kappa],clrs(2,:))
patch([kappa_no_nan,fliplr(kappa_no_nan)],[ymax+0*delte_no_nan,fliplr(delte_no_nan)],clrs(1,:))
patch(kappa([1 idxc_inactive idxc_inactive 1]),repelem(ylim,2),clrs(3,:))

idxc_inactive_theo = find(isfinite(dc_theory),1);
loglog(kappa,dc_theory,'w','LineWidth',2);
loglog(kappa([idxc_inactive_theo idxc_inactive_theo]),[dc_theory(idxc_inactive_theo) ymax ],'w','LineWidth',2);

scatter(kappa,DeltaCMat(2,:),'ow','filled');

kc_inactive = 2/w-1;
loglog(kc_inactive*[1 1],ylim,'w','LineWidth',2)


