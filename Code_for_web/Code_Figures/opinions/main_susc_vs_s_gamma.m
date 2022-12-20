clear; clc; close all
format compact

folder = '../../';
addpath([folder,'functions']);

global A Anw x0 Factor_Xeff List lambda NameOfModel a h
tic

rng(0)

% model of dynamics
NameOfModel = 'opinions'; %  'Neural'; % 'Eco'; %   'Ising-Sch' ; % 'Simple'; %  'SIS';% 'MAK';%  'PD';%


clrs = autumn(10);
clr = clrs(2,:);
mrkrs = 'osd';

gamma = 1;

Sc = 0.99; % 0.972; % 0.9892; %
M = KindOfDynamics( NameOfModel );
range_sol = [0,10]; % range of solutions for high and low states in x
conditions.type = 'BC';
conditions.free_value = 1;
Delta = 10;
N = 1e4;

ChooseExcNodes = 'rand'; % options - {'low degree', 'high degree', 'local', 'rand'}
NetStruct = 'ER';
k0Vec = 9; %round(logspace(log10(10),log10(20),4));
krealVec = zeros(size(k0Vec));
k2 = krealVec; % second moment
all_or_eff = 'eff_F'; % what xss we need: 'all'/ 'eff' /'eff_F'
rhoVec = [0 1e-4];
releaseVec = [0]; % to release or not after holding
reals = 1; % number of realizations

for idxk = 1:length(k0Vec)
    k0 = k0Vec(idxk)

    for idx_release = 1:length(releaseVec)
        release = releaseVec(idx_release);

        for irho = 1: length(rhoVec)
            rho = rhoVec(irho)
            NumExcited = round(rho*N);

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

                toplot = 0; % 'time';

                kappa = k2(idxk)/krealVec(idxk);
                omegaVec = (Sc+logspace(-2,-1/2,7))/kappa;
                Svec = omegaVec*kappa;
                for iw = 1: length(omegaVec)
                    iw
                    omega = omegaVec(iw);
                    S = omega*kappa;
                    xFreals(ir,iw) = find_x_via_omega(omega,M,all_or_eff,toplot,conditions,release);
                end
            end
            xF(irho,:) = mean(xFreals,1);
        end
    end
end

toc

chi = (xF(2,:) - xF(1,:))/(rhoVec(2)-rhoVec(1)) ;



%% Mean Field theory and figures

figure; hold on;

% Sc = 0.9892;

% x = logspace(-6,0,1e3);
% S_th = zeros(2,length(x));
%
% for i=1:2
%     rho = rhoVec(i);
%     S_th(i,:) = -M{1}(x)./ M{2}(x) ./ ((1-rho)*M{3}(x)+rho*M{3}(Delta));
%     x_th(i,:) = interp1(S_th(i,:),x,Svec);
% end
% x_th(isnan(x_th)) = 0;
% chi_th = (x_th(2,:) - x_th(1,:))/(rhoVec(2)-rhoVec(1)) ;



scatter(Svec-0.98,chi,100,clr,'Marker',mrkrs(1),'LineWidth',2);
% plot(S-1,x,'-','Color',colors(1,:),'LineWidth',2);
%     scatter(Svec-1,chi_th(ia,:),100,clrs(ia,:));
plot(Svec-Sc,1.2*(Svec-Sc).^-gamma,'LineWidth',2,'Color',clr)
% title(['$a=',num2str(a,2),'\ \gamma=2-1/a=',num2str(2-1/a,2),'$'],'Interpreter','latex','FontWeight','normal')

axis square tight;
set(gca,'FontSize',20,'linewidth',2,'box','on','XScale','log','YScale','log')
xticks(10.^[-3 -2 -1]); yticks(10.^(-3:3));
axis(10.^[-2 -0.4 0.5 2.3])
ylabel('\boldmath$\chi(\mathcal{S},\rho=0)$','Interpreter','latex','FontSize',25);
xlabel('\boldmath$\mathcal{S}-\mathcal{S}_c$','Interpreter','latex','FontSize',25);



%% save
filename = ['figs/gamma_x_vs_S_','kappa',num2str(round(kappa)),'_N',num2str(N)];
saveas(gcf,filename)

