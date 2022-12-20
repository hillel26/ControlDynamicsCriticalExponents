clear; clc; % close all
format compact

folder = '../../';
addpath([folder,'functions']);

global A Anw x0 Factor_Xeff List lambda NameOfModel a h
tic

% model of dynamics
NameOfModel = 'MM'; %  'Neural'; % 'Eco'; %   'Ising-Sch' ; % 'Simple'; %  'SIS';% 'MAK';%  'PD';%
a = 1;
h = a; % parameters of MM
M = KindOfDynamics( NameOfModel );
range_sol = [0,10]; % range of solutions for high and low states in x
conditions.type = 'BC';
conditions.free_value = 1;
Delta = 10;
N=1e4;

diff_conds = {'free','forced'};
ChooseExcNodes = 'rand'; % options - {'low degree', 'high degree', 'local', 'rand'}
all_or_eff = 'eff_F'; % what xss we need: 'all'/ 'eff' /'eff_F'
rho = 0.03;
NExcVec = round(N*[0,rho]); % round(N*logspace(-4,-1,30)); %
fracVec = NExcVec/N;
releaseVec = [0]; % to release or not after holding

NetStruct = 'ER';
k0Vec = 9; %round(logspace(log10(10),log10(20),4));
krealVec = zeros(size(k0Vec));
k2 = krealVec; % second moment


reals = 1; % number of realizations

colors = [0.6 0 0; 0 0.3 0.6]; % lines(2);
mrkrs = 'osv';

for idxk = 1:length(k0Vec)
    k0 = k0Vec(idxk)
    
    for idx_release = 1:length(releaseVec)
        release = releaseVec(idx_release);
        
        for idx_conds = 1:length(diff_conds)
            conds = diff_conds{idx_conds}
            switch conds
                case 'free'
                    NumExcited = NExcVec(1);
                case 'forced'
                    NumExcited = NExcVec(2);
            end
            
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
                wVec = linspace(0.5,2,8)/kappa;
                for iw = 1:length(wVec)
                    iw
                    omega = wVec(iw);
                    if a>=1
                        toplot = 0;
                        xFreals(ir,iw) = find_x_via_omega(omega,M,all_or_eff,toplot,conditions,release);
                    else
                        
                        A = omega*Anw;
                        
                        x0(x_free) = conditions.free_value;
                        
                        % solve the odes
                        
                        M0 = M{1}; M1 = M{2}; M2 = M{3};
                        
                        dxdt = @(t,x) (M0(x) + M1(x).*(A*M2(x))).*x_free;
                        
                        Tmax = 3e2; % 5e4;
                        dt = 1e-2;
                        t = 0;
                        xt = x0;
                        x = xt;
                        delta = 1; % in order to run first iteration
                        i = 0;
                        while   t(end) <= Tmax && delta > 1e-5
                            i = i+1;
                            % Runga kutta
                            xth = xt + dt/2*dxdt(i*dt,xt);
                            xth(xth<0) = 0;
                            xt = xt + dt*dxdt((i+1/2)*dt,xth);
                            xt(xt<0) = 0;
                            if not(mod(i,50))
                                x = [x,xt];
                                t = [t,i*dt];
                                [delta, ~] = max( abs( (x(end,:) - x(end-1,:)) / ( t(end)-t(end-1) ) ) );
                            end
                        end
                        
%                         figure; plot(t,x)
                        Factor_Xeff_F = sum(A(x_free,x_free),1)/sum(sum(A(x_free,x_free)));
                        xFreals(ir,iw) = Factor_Xeff_F * xt(x_free); % Xeff
                    end
                end
                SReals(ir,:) = wVec*kappa;
            end
            xF(idx_conds,:) = mean(xFreals,1);
            SVec(idx_conds,:) = mean(SReals,1);
        end
    end    
end

toc


%% Mean Field for free system with phases colors
% for given dynamics
figure; hold on;

x = logspace(-5,3,1e3);
i=0;
for e=[0,rho]
    i=i+1;
    beta_func = @(x) -M{1}(x)./ M{2}(x)./((1-e)*M{3}(x)+e*M{3}(Delta));
    beta = beta_func(x);
    p(i) = scatter(SVec(i,:),xF(i,:),100,colors(i,:),'Marker',mrkrs(i),'linewidth',2);
    plot([0 beta],[0 x],'-','Color',colors(i,:),'LineWidth',2);
end

legend(p,{'$\rho=0$',['$\rho=',num2str(rho,2),'$']},'Interpreter','latex','Box','off','Location','northwest')

axis ([0.3, 2.2, -0.2, 2]); axis square;
set(gca,'FontSize',20,'linewidth',2,'box','on','XScale','lin','YScale','lin')
xticks(0:6); yticks(0:6);
ylabel('\boldmath$\bar{\rm x}$','Interpreter','latex','FontSize',25);
xlabel('\boldmath$\mathcal{S}$','Interpreter','latex','FontSize',25);


%% save
filename = ['figs/x_vs_S_',NetStruct,'_kappa',num2str(round(kappa)),'_N',num2str(N)];
saveas(gcf,filename)

