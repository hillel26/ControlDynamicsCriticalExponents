clear; clc; close all;
format compact

global a h
addpath('../../functions');

NameOfModel = 'MM';
M = KindOfDynamics(NameOfModel);

%% colors

bi_color = 0.3*[1 1 1];
inactive_color = [0.6 0 0]; %[155 0 0]/255;
active_color = [0 0.3 0.6];
forced_color = [0 0.6 0];

%% MM illust free and forced
a = 2; h = 2;
Delta = 10;
x = linspace(0,12,801);

% free
figure; hold on

rho = 0;   
beta_free = x.^a./((1-rho)./(1+x.^-h) + rho./(1+Delta^-h));

plot(beta_free,x,'Color',active_color,'LineWidth',5)
plot([0 1.025],[0 0],'Color',inactive_color,'LineWidth',5)

set(gca,'XTick',[],'YTick',[],'LineWidth',2,'Layer','top','FontSize',20)
axis([0,2.5,-0.3,1.5]); axis square; box on


% forced
figure; hold on

rho = 2e-3;   
beta_forced = x.^a./((1-rho)./(1+x.^-h) + rho./(1+Delta^-h));

plot(beta_free,x,'Color',[active_color 0.3],'LineWidth',5)
plot([0 1.025],[0 0],'Color',[inactive_color 0.3],'LineWidth',5)
plot(beta_forced,x,'Color',forced_color,'LineWidth',5)

set(gca,'XTick',[],'YTick',[],'LineWidth',2,'Layer','top','FontSize',20)
axis([0,2.5,-0.3,1.5]); axis square; box on


%% cases of a and h

a = 1; h = 2;
x = linspace(0,12,801);

% free
figure; hold on

beta_free = x.^a.*(1+x.^-h);

plot(beta_free(x>=1),x(x>=1),'Color',active_color,'LineWidth',5)
plot(beta_free(x<1),x(x<1),'--','Color',active_color,'LineWidth',5)
plot([0 5],[0 0],'Color',inactive_color,'LineWidth',5)

set(gca,'XTick',[],'YTick',[],'LineWidth',2,'Layer','top','FontSize',20)
axis([0,5,-0.3,3]); axis square; box on


a = 2; h = 1;
x = linspace(0,12,801);

% free
figure; hold on

rho = 0;   
beta_free = x.^a.*(1+x.^-h) ;

plot(beta_free,x,'Color',active_color,'LineWidth',5)
plot([0 5],[0 0],'--','Color',inactive_color,'LineWidth',5)

set(gca,'XTick',[],'YTick',[],'LineWidth',2,'Layer','top','FontSize',20)
axis([0,5,-0.3,3]); axis square; box on
