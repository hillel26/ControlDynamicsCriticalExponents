clear; format compact;
% close all; clc;
fold = 'c:\Users\user\Dropbox\Havlin\DynamicsRecovering\single_node_recovery\code\';
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
wVec = logspace(-1,2,20);
k0 = 4;
lmax = 3; % only for tree
kRealVec = zeros(size(wVec));
k2 = kRealVec; % second moment

free_value_vec = [0 10];
NExcVec = [0 1]; %unique( round(logspace(0,log10(NVec/25),21)) );
DeltaCMat = zeros(1,length(wVec));
DeltaCMatstd = DeltaCMat;
Critical_delta = 1;
reals = 5;
CrDelRealVec = zeros(1,reals);

for idxN = 1:length(NVec)
    N = NVec(idxN)
    for idxw = 1:length(wVec)
        w = wVec(idxw)
        % build the network
        for idxl = 1:length(lambdaVec)
            NetStruct = 'ER';
            parameters = k0;
            if strcmp(NetStruct,'SF')
                lambda = lambdaVec(idxl);
                parameters = [lambda w];
            end
            
            %     if strcmp(ChooseExcNodes,'local')
            %         List = MatrixToList(Anw);
            %     end
            
            % initial/boundary condition
            
            holding_value = 1;
            
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
                    if strcmp(NetStruct,'RT'); N=1+w*((w-1)^lmax - 1)/(w-2); end
                    Anw = BuildNetwork(N, NetStruct, parameters,'gcc'); % Adjacency matrix, not weighted
                    A = w*Anw;
                    Factor_Xeff = sum(A,1)/sum(sum(A));

                    n = size(Anw,1);
                    if strcmp(NetStruct,'SF')
                        kRealVec(idxl) = mean(sum(Anw)); k2(idxl) = mean(sum(Anw).^2);
                    else
                        kRealVec(idxw) = mean(sum(Anw));
                        k2(idxw) = mean(sum(Anw).^2);
                    end
                    x0 = Initial_or_fixed_Conditions(n,NumExcited,ChooseExcNodes,conditions,holding_value);
                    toplot = 0;
                    
                    % find critical omega
                    options = optimset('TolX',1e-3);
                    %                         check limits
                    xA = find_x_via_delta(1e-3,M,all_or_eff,toplot,conditions,1);
                    xB = find_x_via_delta(20,M,all_or_eff,toplot,conditions,1);
                    if xB<x_th
                        Critical_delta = inf;
                    elseif xA>x_th
                        Critical_delta = nan;
                    else
                        Critical_delta = ...
                            fzero(@(delta) find_x_via_delta(delta,M,all_or_eff,toplot,conditions,1) - x_th ,[1e-3,20] ,options);
                    end
                    
                    % PlotExamples when ...
                    for Examples=[]
                        if  (idxNE==1) && (ir==1) && (idxl==1) && (idxN==length(NVec)) && (idxw==1)
                            AxXvec = linspace(eps,2*Critical_delta,50);
                            AxYvec = arrayfun(@(omega) find_x_via_omega(omega,x0,M,all_or_eff,toplot,conditions,release),AxXvec);
                            arrayfun(@(omega) find_x_via_omega(omega,x0,M,'eff',1,conditions,release),Critical_delta*[0.99,1.01]);
                            
                            figure; % Xeff vs omega
                            plot(AxXvec, AxYvec,'.'...
                                ,AxXvec,(Factor_Xeff*x0)*ones(size(AxXvec)),'.')
                            xlabel('\omega_c'); ylabel('x_{eff}'); legend('SteadyState','InitialCondition');
                            title(['k = ',num2str(kRealVec(idxw)),', \rho = ',num2str(num2str(NExcVec(idxNE)/n))]);
                            set(gca,'XTick',sort([xlim,Critical_delta]),'XGrid','on','box','on','fontsize',15);
                            
                            figure; % Xeff vs beta
                            plot(AxXvec*k2(idxw)/kRealVec(idxw),AxYvec,'.'...
                                ,AxXvec*k2(idxw)/kRealVec(idxw),(Factor_Xeff*x0)*ones(size(AxXvec)),'.')
                            xlabel('\beta_c'); ylabel('x_{eff}'); legend('SteadyState','InitialCondition');
                            title(['k = ',num2str(kRealVec(idxw)),', \rho = ',num2str(num2str(NExcVec(idxNE)/n))]);
                            set(gca,'fontsize',15); box on
                        end
                    end
                    
                    % X0eff = Factor_Xeff*x0;
                    % beta = omega*sum(sum(A^2))/sum(A(:));
                    CrDelRealVec(ir) = Critical_delta;
                end
                DeltaCMat(idx_cond, idxw) = mean(CrDelRealVec);
                DeltaCMatstd(idx_cond,idxw) = std(CrDelRealVec);
                
            end
        end
    end
    toc
end

kappa = k2./kRealVec-1;
kappa = mean(kappa);

% theory of intersection between functions
% w_theory = logspace(-1/2,2,21);
dc_theory = arrayfun(@(w) find_dc_by_wk_theory_intersection(w,kappa,range) ,  wVec);

%% Figure  vs k 
figure; hold on

set(gca,'FontSize',20,'box','on','LineWidth',1.5,'Layer','top','XScale','log','YScale','log')
ylabel('\boldmath$\Delta$','Interpreter','latex','FontSize',33)
xlabel('\boldmath$w$','Interpreter','latex','FontSize',33)
axis square

ymin = 1e-3; ymax = 5;
xlim(wVec([1,end]))
ylim([ymin,ymax])

% patches
clrs = [0 51 102; 255 192 0; 155 0 0]/255;
w_no_nan = wVec(isfinite(DeltaCMat(2,:)));
delta_no_nan = DeltaCMat(2,isfinite(DeltaCMat(2,:)));
idxc_rec = find(isnan(DeltaCMat(1,:)),1);

patch([wVec,fliplr(wVec)],[ymin+0*wVec,ymax+0*wVec],clrs(2,:))
patch([w_no_nan,fliplr(w_no_nan)],[ymax+0*delta_no_nan,fliplr(delta_no_nan)],clrs(1,:))
patch(wVec([1, idxc_rec, idxc_rec, 1]),repelem(ylim,2), clrs(3,:))

idxc_rec_theo = find(isfinite(dc_theory),1);
loglog(wVec,dc_theory,'w','LineWidth',2);
loglog(wVec(idxc_rec_theo)*[1 1],[dc_theory(idxc_rec_theo) ymax],'w','LineWidth',2);

scatter(wVec,DeltaCMat(2,:),'ow','filled');

wc_inactive = 2/(kappa+1);
loglog(wc_inactive*[1 1],ylim,'w','LineWidth',2);

