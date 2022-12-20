clear; clc;  close all
folder = '../../';
addpath([folder,'functions']);

global A Anw x0 Factor_Xeff List lambda NameOfModel a h
tic

% model of dynamics
NameOfModel = 'SIS'; %  'Neural'; % 'Eco'; %   'Ising-Sch' ; % 'Simple'; %  'SIS';% 'MAK';%  'PD';%

clr = summer(10);
clr = clr(2,:);
mrkrs = 'osd';

delta = 2;

M = KindOfDynamics( NameOfModel );
range_sol = [0,10]; % range of solutions for high and low states in x
conditions.type = 'BC';
conditions.free_value = 1;
Delta = 10;
N=1e4;
Sc = 0.99;

ChooseExcNodes = 'rand'; % options - {'low degree', 'high degree', 'local', 'rand'}
NetStruct = 'ER';
k0Vec = 9; %round(logspace(log10(10),log10(20),4));
krealVec = zeros(size(k0Vec));
k2 = krealVec; % second moment
all_or_eff = 'eff_F'; % what xss we need: 'all'/ 'eff' /'eff_F'

rhoVec = logspace(-3,-1,10);
releaseVec = [0]; % to release or not after holding
reals = 1; % number of realizations

for idxk = 1:length(k0Vec)
    k0 = k0Vec(idxk)

    for idxrho = 1:length(rhoVec)
        rho = rhoVec(idxrho)
        NumExcited = round(N*rho); % round(N*logspace(-4,-1,30)); %

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
                x_free = not(logical(x0));

                kappa = k2(idxk)/krealVec(idxk);
                omega = Sc/kappa;

                toplot = 0; %'time';
                xFreals(ir) = find_x_via_omega(omega,M,all_or_eff,toplot,conditions,release);

                SReals(ir) = omega*kappa;
            end
            xF(idxrho) = mean(xFreals,1);
            SVec(idxrho) = mean(SReals,1);
        end
    end
end

toc


%% Mean Field for free system with phases colors
% for given dynamics


% x_theory = 0*rhoVec;
% i=0;
% for rho = rhoVec
%     i=i+1;
%     beta_func = @(x) -M{1}(x)./ M{2}(x)./((1-rho)*M{3}(x) + (rho*M{3}(Delta)));
%     x_theory(i) = fzero(@(x) beta_func(x)-1,rho^(1/delta));
% end

%%

figure; hold on

scatter(rhoVec,xF,100,clr,'LineWidth',2,'Marker',mrkrs(1));
% loglog(rhoVec,x_theory,'-','LineWidth',2);
loglog(rhoVec,1.6*rhoVec.^(1/delta),'-','LineWidth',2,'Color',clr)

axis square tight;
axis(10.^[-3 -1 -1.5 0])

set(gca,'FontSize',20,'box','off','linewidth',2,'box','on','XScale','log','YScale','log')
xticks(10.^[-3 -2 -1]); yticks(10.^[-3 -2 -1 0]);
ylabel('\boldmath$\bar{\rm x}(\mathcal{S}_c,\rho)$','Interpreter','latex','FontSize',25);
xlabel('\boldmath$\rho$','Interpreter','latex','FontSize',25);



%% save
filename = ['delta_x_vs_rho_',NetStruct,'_kappa',num2str(round(kappa)),'_N',num2str(N)];
saveas(gcf,filename)

