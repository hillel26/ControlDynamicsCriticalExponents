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

%% SIS illust free and forced

NameOfModel = 'SIS'; %  'Neural'; % 'Eco'; %   'Ising-Sch' ; % 'Simple'; %  'SIS';% 'MAK';%  'PD';%
M = KindOfDynamics( NameOfModel );

Delta = 10;
x = linspace(0,1,801);

% free
figure; hold on

rho = 0;  
beta_free =  -M{1}(x)./ M{2}(x)./((1-rho)*M{3}(x)+rho*M{3}(Delta));

plot(beta_free,x,'Color',active_color,'LineWidth',5)
plot([0 1.025],[0 0],'Color',inactive_color,'LineWidth',5)

set(gca,'LineWidth',2,'Layer','top','FontSize',20)
axis([0,2.5,-0.3,1.5]); axis square; box on
xticks([]); yticks([])

% forced
figure; hold on

rho = 2e-3;   
beta_forced =  -M{1}(x)./ M{2}(x)./((1-rho)*M{3}(x)+rho*M{3}(Delta));

plot(beta_free,x,'Color',[active_color 0.3],'LineWidth',5)
plot([0 1.025],[0 0],'Color',[inactive_color 0.3],'LineWidth',5)
plot(beta_forced,x,'Color',forced_color,'LineWidth',5)

set(gca,'XTick',[],'YTick',[],'LineWidth',2,'Layer','top','FontSize',20)
axis([0,2.5,-0.3,1.5]); axis square; box on



%% SIS illust free supplementary

NameOfModel = 'SIS'; %  'Neural'; % 'Eco'; %   'Ising-Sch' ; % 'Simple'; %  'SIS';% 'MAK';%  'PD';%
M = KindOfDynamics( NameOfModel );

x = linspace(-1,1,801);

% free
figure; hold on

beta_free =  -M{1}(x)./ M{2}(x)./M{3}(x);

plot(beta_free(x>=0),x(x>=0),'Color',active_color,'LineWidth',5)
plot(beta_free(x<-0.05),x(x<-0.05),'--','Color',active_color,'LineWidth',5)

plot([0 1.025],[0 0],'Color',inactive_color,'LineWidth',5)
plot([1.1 5],[0 0],'--','Color',inactive_color,'LineWidth',5)


set(gca,'LineWidth',2,'Layer','top','FontSize',20)
axis([0,5,-0.5,1.5]); axis square; box on
xticks([]); yticks([])

