clear; clc; close all
format compact
folder = '../../';
addpath([folder,'functions']);

global A Anw x0 Factor_Xeff List lambda NameOfModel a h
tic


rng(0)

% model of dynamics
NameOfModel = 'MM'; %  'Neural'; % 'Eco'; %   'Ising-Sch' ; % 'Simple'; %  'SIS';% 'MAK';%  'PD';%

leginfo = {};
aVec = [1/2 1 2];
clrs = parula(10);
clrs = clrs(2:4,:);
mrkrs = 'osd';
for ia = 1:length(aVec)
a = aVec(ia);
h=a; % parameters of MM
beta = 1/a;
Sc = 0.9892;
M = KindOfDynamics( NameOfModel );
range_sol = [0,10]; % range of solutions for high and low states in x
conditions.type = 'BC';
Delta = 10;
N=1e3;

diff_conds = {'free low'}; %{'free low','free high','ignited'};

ChooseExcNodes = 'rand'; % options - {'low degree', 'high degree', 'local', 'rand'}
NetStruct = 'ER';
k0Vec = 9; %round(logspace(log10(10),log10(20),4));
krealVec = zeros(size(k0Vec));
k2 = krealVec; % second moment
all_or_eff = 'eff_F'; % what xss we need: 'all'/ 'eff' /'eff_F'
free_value_vec = [1e-2, 1]; % the initial value of the free nodes
rho = 0.05;
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
            omegaVec = (Sc+logspace(-3,-1/2,15))/kappa;
            Svec = omegaVec*kappa;
            for iw = 1: length(omegaVec)
                iw
                omega = omegaVec(iw);
                S = omega*kappa;
                conditions.free_value = (S-Sc).^beta;
                xFreals(ir,iw) = find_x_via_omega(omega,M,all_or_eff,toplot,conditions,release);
            end
            SReals(ir,:) = omegaVec*k2(idxk)/krealVec(idxk);
        end
        xF(ia,:) = mean(xFreals,1);
        SVec = mean(SReals,1);
    end
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

for ia=1:length(aVec)
a = aVec(ia);
beta = 1/a;

p(ia) = scatter(SVec-Sc,xF(ia,:),100,clrs(ia,:),'linewidth',2,'Marker',mrkrs(ia));
plot(SVec-1,0.7*(SVec-1).^beta,'LineWidth',2,'Color',clrs(ia,:))

leginfo{ia} = ['\beta=',num2str(beta)];
end



axis square tight;
set(gca,'FontSize',20,'linewidth',2,'box','on','XScale','log','YScale','log')
xticks(10.^[-3 -2 -1]); yticks(10.^(-6:2:0));
% axis(10.^[-3 -1 -6.5 -2])
ylabel('\boldmath$\bar{\rm x}(\mathcal{S},\rho=0)$','Interpreter','latex','FontSize',25);
xlabel('\boldmath$\mathcal{S}-\mathcal{S}_c$','Interpreter','latex','FontSize',25);
% title(['$a=',num2str(a,2),'\ \beta=1/a=',num2str(1/a,2),'$'],'Interpreter','latex','FontWeight','normal')
legend(p,leginfo,'Box','off','Location','southeast')


%% save
filename = ['figs/beta_x_vs_S_','kappa',num2str(round(kappa)),'_N',num2str(N)];
saveas(gcf,filename)

