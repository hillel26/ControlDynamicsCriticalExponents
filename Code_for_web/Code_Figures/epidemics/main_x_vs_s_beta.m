clear; clc; close all
format compact
folder = '../../';
addpath([folder,'functions']);

global A Anw x0 Factor_Xeff List lambda NameOfModel a h
tic


rng(0)

% model of dynamics
NameOfModel = 'SIS'; %  'Neural'; % 'Eco'; %   'Ising-Sch' ; % 'Simple'; %  'SIS';% 'MAK';%  'PD';%

clr = summer(10);
clr = clr(2,:);
mrkrs = 'osd';

beta = 1;
Sc = 0.99;
M = KindOfDynamics( NameOfModel );
range_sol = [0,10]; % range of solutions for high and low states in x
conditions.type = 'BC';
conditions.free_value = 1/2;
Delta = 10;
N=1e4;

ChooseExcNodes = 'rand'; % options - {'low degree', 'high degree', 'local', 'rand'}
NetStruct = 'ER';
k0Vec = 9; %round(logspace(log10(10),log10(20),4));
krealVec = zeros(size(k0Vec));
k2 = krealVec; % second moment
all_or_eff = 'eff_F'; % what xss we need: 'all'/ 'eff' /'eff_F'

NumExcited = 0;
releaseVec = [0]; % to release or not after holding
reals = 1; % number of realizations


for idxk = 1:length(k0Vec)
    k0 = k0Vec(idxk)
    
    for idx_release = 1:length(releaseVec)
        release = releaseVec(idx_release);
                
        for ir = 1:reals
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
            x0 = Initial_or_fixed_Conditions(n,NumExcited,ChooseExcNodes,conditions,Delta);
            toplot = 0; % 'time';
            
            kappa = k2(idxk)/krealVec(idxk);
            omegaVec = (Sc+logspace(-2,-1/2,10))/kappa;
            Svec = omegaVec*kappa;
            for iw = 1: length(omegaVec)
                iw
                omega = omegaVec(iw);
                S = omega*kappa;
                xFreals(ir,iw) = find_x_via_omega(omega,M,all_or_eff,toplot,conditions,release);
            end
            SReals(ir,:) = omegaVec*k2(idxk)/krealVec(idxk);
        end
        xF = mean(xFreals,1);
        SVec = mean(SReals,1);
    end
end

toc


%% plot with scaling and theory

figure; hold on;

% Sc = 0.9892;

% x = logspace(-6,3,1e3);
% S_func = @(x) -M{1}(x)./ M{2}(x) ./ M{3}(x);
% S = S_func(x);
% plot(S-1,x,'-','Color',clrs(ia,:),'LineWidth',2);


scatter(SVec-0.99,xF,100,clr,'linewidth',2,'Marker',mrkrs(1));
plot(SVec-1,0.7*(SVec-1).^beta,'LineWidth',2,'Color',clr)


axis square tight;
set(gca,'FontSize',20,'linewidth',2,'box','on','XScale','log','YScale','log')
xticks(10.^[-3 -2 -1]); yticks(10.^(-6:1:0));
axis(10.^[-2 -0.5 -2.3 -0.5])
ylabel('\boldmath$\bar{\rm x}(\mathcal{S},\rho=0)$','Interpreter','latex','FontSize',25);
xlabel('\boldmath$\mathcal{S}-\mathcal{S}_c$','Interpreter','latex','FontSize',25);
% title(['$a=',num2str(a,2),'\ \beta=1/a=',num2str(1/a,2),'$'],'Interpreter','latex','FontWeight','normal')


%% save
filename = ['figs/beta_x_vs_S_','kappa',num2str(round(kappa)),'_N',num2str(N)];
saveas(gcf,filename)

