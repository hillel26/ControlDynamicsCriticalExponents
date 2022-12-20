clear; format compact;
% close all; clc;
fold = 'c:\Users\user\Dropbox\Havlin\DynamicsRecovering\single_node_recovery\code\';
addpath([fold,'functions']);

global A Anw x0 Factor_Xeff List lambda NameOfModel a b

tic

% model of dynamics
NameOfModel = 'MM'; % 'Eco'; % 'Simple'; % 'Glauber'; % 'Voter'; % 'SIS';% 'MAK';%  'PD';%
a = 1; % exponents of MM dynamics

h = gobjects(0);

M = KindOfDynamics( NameOfModel );
range = [0,10]; % range of values of x for theory
ChooseExcNodes = 'rand'; % options - {'low degree', 'high degree', 'local', 'rand','1'}
conditions.type = 'BC';% 'IC'; %

conds = {'free high','ignite'};
release = 1;
all_or_eff = 'eff'; % what xss we need: 'all'/ 'eff' /'eff_F'/ 'average'

x_th = 1/2;
HoldingValueVec = 1; % logspace(-1,2,4);
NVec = 5e3; %round(logspace(1,4,11)); %1e5*(1:10);
lambdaVec = 3; % linspace(2.2,3.2,5);
k0Vec = logspace(log10(3),log10(300),15);
lmax = 3; % only for tree
kRealVec = zeros(size(k0Vec));
k2 = kRealVec; % second moment
free_value_vec = [0 10];
NExcVec = [0 0.01*NVec]; %unique( round(logspace(0,log10(NVec/25),21)) );
OmegaCMat = zeros(length(HoldingValueVec),length(k0Vec));
OmegaCMatstd = OmegaCMat;
Critical_omega = 1;
reals = 10;
CrOmRealVec = zeros(1,reals);

for idxN = 1:length(NVec)
    N = NVec(idxN)
    for idxk = 1:length(k0Vec)
        k0 = k0Vec(idxk)
        % build the network
        for idxl = 1:length(lambdaVec)
            NetStruct = 'ER';
            parameters = k0;
            if strcmp(NetStruct,'SF')
                lambda = lambdaVec(idxl)
                parameters = [lambda k0];
            end
            
            %     if strcmp(ChooseExcNodes,'local')
            %         List = MatrixToList(Anw);
            %     end
            
            % initial/boundary condition
            
            for idxval = 1:length(HoldingValueVec)
                holding_value = HoldingValueVec(min(idxval,end))
                
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
                        %                         check limits
                        %                         find_x_via_omega(1e-2,M,all_or_eff,toplot,conditions,1)
                        %                         find_x_via_omega(20,M,all_or_eff,toplot,conditions,1)
                        %                         x_th
                        Critical_omega = ...
                            fzero(@(omega) find_x_via_omega(omega,M,all_or_eff,toplot,conditions,1) - x_th ,[1e-3,100] ,options);
                        
                        % PlotExamples when ...
                        for Examples=[]
                            if  (idxNE==1) && (ir==1) && (idxl==1) && (idxN==length(NVec)) && (idxk==1)
                                AxXvec = linspace(eps,2*Critical_omega,50);
                                AxYvec = arrayfun(@(omega) find_x_via_omega(omega,x0,M,all_or_eff,toplot,conditions,release),AxXvec);
                                arrayfun(@(omega) find_x_via_omega(omega,x0,M,'eff',1,conditions,release),Critical_omega*[0.99,1.01]);
                                
                                figure; % Xeff vs omega
                                plot(AxXvec, AxYvec,'.'...
                                    ,AxXvec,(Factor_Xeff*x0)*ones(size(AxXvec)),'.')
                                xlabel('\omega_c'); ylabel('x_{eff}'); legend('SteadyState','InitialCondition');
                                title(['k = ',num2str(kRealVec(idxk)),', \rho = ',num2str(num2str(NExcVec(idxNE)/n))]);
                                set(gca,'XTick',sort([xlim,Critical_omega]),'XGrid','on','box','on','fontsize',15);
                                
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
                        CrOmRealVec(ir) = Critical_omega;
                    end
                    if strcmp(NetStruct,'SF')
                        OmegaCMat(idx_cond,idxl) = mean(CrOmRealVec);
                        OmegaCMatstd(idx_cond,idxl) = std(CrOmRealVec);
                    else
                        OmegaCMat(idx_cond,idxk) = mean(CrOmRealVec);
                        OmegaCMatstd(idx_cond,idxk) = std(CrOmRealVec);
                    end
                end
            end
        end
    end
    toc
end

BetaCMat = OmegaCMat*diag(k2./kRealVec);
kappa = k2./kRealVec-1;

%% theory ignite fraction

e = NExcVec(2)/NVec;
Delta = HoldingValueVec;
x = linspace(range(1),range(2),1e3);
beta_theory = -M{1}(x)./M{2}(x)./(M{3}(x)+e*M{3}(Delta));
[bc_theory,xc] = findpeaks(beta_theory,x);
wc_theory = bc_theory./(kappa+1);


%% theory ignite single
options = optimset('TolX',1e-3);
wc_theory_single = arrayfun(@(kappa) fzero(@(w) is_rcoverable_theory(w,kappa,holding_value,range) - 1/2 , [2/(kappa+1),1e1],options) , kappa);


%%
figure; hold on
set(gca,'FontSize',20,'box','on','LineWidth',2,'XScale','log','YScale','log','layer','top')
xlabel('\boldmath$\kappa$','Interpreter','latex','FontSize',33)
ylabel('\boldmath$w$','Interpreter','latex','FontSize',33)
axis square

ymin = 10^-2.5; ymax = 4;
xlim([min(kappa),max(kappa)])
ylim([ymin,ymax])

% patches
clrs = [100 170 255; 255 192 0; 155 0 0]/255; % [0 51 102; 255 192 0; 155 0 0]/255;

patch([xlim,fliplr(xlim)],repelem(ylim,2),clrs(2,:))
patch([kappa,fliplr(kappa)],[ymax+0*OmegaCMat(2,:),fliplr(OmegaCMat(2,:))],clrs(1,:))
patch([kappa,fliplr(kappa)],[ymin+0*OmegaCMat(1,:),fliplr(OmegaCMat(1,:))],clrs(3,:))
patch([kappa,fliplr(kappa)],[ymax+0*wc_theory_single,fliplr(wc_theory_single)],[0 51 102]/255)


% inactive phase
wc_inactive_theory = 2./(kappa+1);

% line and symbols
loglog(kappa,wc_theory,'w','LineWidth',2)
% loglog(kappa,wc_theory_single,'w','LineWidth',2)
% loglog(kappa,wc_inactive_theory,'w','LineWidth',2)
% scatter(kappa, OmegaCMat(1,:),'wo','filled')
% scatter(kappa, OmegaCMat(2,:),'wo','filled')
xticks([1 10 100])



