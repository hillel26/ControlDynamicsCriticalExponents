clear; clc; % close all
folder = 'C:/Users/user/Dropbox/Havlin/DynamicsRecovering/single_node_recovery/code/';
addpath([folder,'functions']);

global A Anw x0 Factor_Xeff List lambda NameOfModel a b
tic

% model of dynamics
NameOfModel = 'MM'; %  'Neural'; % 'Eco'; %   'Ising-Sch' ; % 'Simple'; %  'SIS';% 'MAK';%  'PD';%
% a=1; b=2; % parameters of MM
M = KindOfDynamics( NameOfModel );
range_sol = [0,10]; % range of solutions for high and low states in x
conditions.type = 'BC';
holding_value = 10;
N=1e4;

diff_conds = {'free low','free high','ignited'};

ChooseExcNodes = 'rand'; % options - {'low degree', 'high degree', 'local', 'rand'}
NetStruct = 'PPI_Yeast';
k0Vec = 2; %round(logspace(log10(10),log10(20),4));
krealVec = zeros(size(k0Vec));
k2 = krealVec; % second moment
all_or_eff = 'eff_F'; % what xss we need: 'all'/ 'eff' /'eff_F'
free_value_vec = [0, 50]; % the initial value of the free nodes
NExcVec = [0 1]; %round(N*[0,0.1]); % round(N*logspace(-4,-1,30)); %
fracVec = NExcVec/N;
releaseVec = [1]; % to release or not after holding
Critical_w = 0.1;
real = 10; % number of realizations
leg = cell(size(k0Vec));
colors =[1 0 0; 0 0 1; 0.7 0 1]; % lines(2);
figure; hold on;
h = gobjects(0); % returns an empty graphics object array.
mrkrs = 'ovs';

for idxk = 1:length(k0Vec)
    k0 = k0Vec(idxk)
    
    for idx_release = 1:length(releaseVec)
        release = releaseVec(idx_release);
        
        for idx_conds = 1:length(diff_conds)
            conds = diff_conds{idx_conds}
            switch conds
                case 'free low'
                    conditions.free_value = free_value_vec(1);
                    NumExcited = NExcVec(1);
                case 'free high'
                    conditions.free_value = free_value_vec(2);
                    NumExcited = NExcVec(1);
                case 'ignited'
                    conditions.free_value = free_value_vec(1);
                    NumExcited = NExcVec(2);
            end
            
            for ir = 1:real
                ir
                % build the network
                %     rng(0);
                parameters = k0;
                if strcmp(NetStruct,'SF')
                    lambda = 3;
                    parameters = [lambda k0];
                end
                Anw = BuildNetwork(N, NetStruct,parameters,'gcc'); % Adjacency matrix, not weighted
                n = size(Anw,1);
                krealVec(idxk) = mean(sum(Anw));
                k2(idxk) = mean(sum(Anw).^2);
                % initial/fixed condition
                x0 = Initial_or_fixed_Conditions(n,NumExcited,ChooseExcNodes,conditions,holding_value);
                toplot=0;
                
                %             % find critical omega
                %             options = optimset('TolX',1e-4);
                %             Critical_omega = ...
                %                 fzero(@(omega) find_x_via_omega(omega,x0,M,all_or_eff,toplot,conditions,1)...
                %                 - 1e-2 ,[eps, 10*Critical_omega] ,options);
                
                %         OmegaVec = linspace(eps,2*Critical_omega,50);
                wVec = logspace(-2,1,100);
                XeffVec = arrayfun(@(omega) find_x_via_omega(omega,M,all_or_eff,toplot,conditions,release),wVec);
                BetaVec = wVec*k2(idxk)/krealVec(idxk);
                %     BetaVec = OmegaVec*sum(sum(Anw^2))/sum(sum(Anw));
            %                 XeffVec = mean(XeffVec);
            %                 BetaVec = mean(BetaVec);
            % Xeff vs beta
            scatter(wVec,XeffVec,50,mrkrs(idx_conds),'filled','markerfacecolor',colors(idx_conds,:));
            plot(wVec,XeffVec,'Color',colors(idx_conds,:));
            end
            %     figure(2)
            %     plot(BetaVec,(Factor_Xeff*x0)*ones(size(BetaVec)),'.')
            
            %                 xeff_vs_rho(idx_rho) = XeffVec;
        end
    end
    
end

kappa = k2(idxk)/krealVec(idxk)-1;

toc

%% Theory for recoverability
options = optimset('TolX',1e-3);
wc_theory = fzero(@(w) is_rcoverable_theory(w,kappa,holding_value,range_sol) - 1/2 , [1,50]/kappa,options);


%% Mean Field for free system with phases colors
% for given dynamics

x = logspace(-2,2,1001);
beta_func = @(x) -M{1}(x)./( M{2}(x).*M{3}(x) );
beta = beta_func(x);
w = beta/kappa;
[bc, xc] = findpeaks(-beta,x); bc = -bc(end);


plot(w,x,'--k','LineWidth',2.2);
plot(w(x>=xc),x(x>=xc),'Color','k','LineWidth',2.2);
plot(wVec,0*wVec,'k','LineWidth',2.2);

plot(wc_theory*[1 1],[0,interp1(w(x>=xc),x(x>=xc),wc_theory)],'--','Color',[0.7 0 1],'LineWidth',2.2);


axis([1e-2, 1e1, 0,10]); axis square;
set(gca,'FontSize',20,'box','off','linewidth',1.5,'box','on','XScale','log','YScale','lin')
ylabel('\boldmath$\bar{x}$','Interpreter','latex','FontSize',30);
xlabel('\boldmath$w$','Interpreter','latex','FontSize',30);

% scatter(BetaVec,XeffVec,'bo','filled');

scs = findobj(gcf,'type','scatter');
legend(flipud(scs([1,11,21])),{'low','high','ignited'},'box','off','Location','northwest')



