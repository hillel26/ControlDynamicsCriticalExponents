clear; clc; close all;
format compact
 
addpath('../../functions');


%% colors

bi_color = 0.3*[1 1 1];
inactive_color = [0.6 0 0]; %[155 0 0]/255;
active_color = [0 0.3 0.6];
forced_color = [0 0.6 0];

%% SIS illust free and forced

NameOfModel = 'opinions'; %  'Neural'; % 'Eco'; %   'Ising-Sch' ; % 'Simple'; %  'SIS';% 'MAK';%  'PD';%
M = KindOfDynamics( NameOfModel );

Delta = 10;
x_free = linspace(-3,3,1e2);
xp = linspace(0,3,1e2);
xn = linspace(-0.45,-0.1,1e2);
xnn = linspace(-3,-0.45,1e2);

% free
figure; hold on

rho = 0;  
beta_free =  -M{1}(x_free)./ M{2}(x_free)./((1-rho)*M{3}(x_free)+rho*M{3}(Delta));

plot(beta_free,x_free,'Color',active_color,'LineWidth',5)
plot([0 1.015],[0 0],'Color',inactive_color,'LineWidth',5)

set(gca,'LineWidth',2,'Layer','top','FontSize',20)
axis([0.5,2,-3,3]); axis square; box on
xticks([]); yticks([])

% forced
figure; hold on

rho = 5e-2;   
beta_forced_p =  -M{1}(xp)./ M{2}(xp)./((1-rho)*M{3}(xp)+rho*M{3}(Delta));
beta_forced_nn =  -M{1}(xnn)./ M{2}(xnn)./((1-rho)*M{3}(xnn)+rho*M{3}(Delta));
beta_forced_n =  -M{1}(xn)./ M{2}(xn)./((1-rho)*M{3}(xn)+rho*M{3}(Delta));

plot(beta_free,x_free,'Color',[active_color 0.3],'LineWidth',5)
plot([0 1.015],[0 0],'Color',[inactive_color 0.3],'LineWidth',5)
plot(beta_forced_p,xp,'Color',forced_color,'LineWidth',5)
plot(beta_forced_n,xn,'--','Color',forced_color,'LineWidth',5)
plot(beta_forced_nn,xnn,'-','Color',forced_color,'LineWidth',5)

set(gca,'LineWidth',2,'Layer','top','FontSize',20)
axis([0.5,2,-3,3]); axis square; box on
% xticks([]); yticks([])


