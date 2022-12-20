clear; clc;  close all
folder = 'C:/Users/user/Dropbox/Havlin/DynamicsRecovering/single_node_recovery/code/';
addpath([folder,'functions']);

global A Anw x0 Factor_Xeff List lambda NameOfModel a b
tic

% model of dynamics
NameOfModel = 'MM'; %  'Neural'; % 'Eco'; %   'Ising-Sch' ; % 'Simple'; %  'SIS';% 'MAK';%  'PD';%
% a=1; b=2; % parameters of MM
M = KindOfDynamics( NameOfModel );
conditions.type = 'BC';
holding_value = 5;
N=1e4;

ChooseExcNodes = 'rand'; % options - {'low degree', 'high degree', 'local', 'rand'}
NetStruct = 'PPI_Human';
k = 2; %round(logspace(log10(10),log10(20),4));
w = 0.2; % weight

all_or_eff = 'eff_F'; % what xss we need: 'all'/ 'eff' /'eff_F'
release = 0; % to release or not after holding
conditions.free_value = 0;
NumExcited = 1;

% build the network
parameters = k;
if strcmp(NetStruct,'SF')
    lambda = 3;
    parameters = [lambda k];
end
Anw = BuildNetwork(N, NetStruct,parameters,'gcc'); % Adjacency matrix, not weighted
n = size(Anw,1);
kreal = mean(sum(Anw));
k2 = mean(sum(Anw).^2);

% initial/fixed condition
x0 = Initial_or_fixed_Conditions(n,NumExcited,ChooseExcNodes,conditions,holding_value);

toplot = 'pics';
find_x_via_omega(w,M,all_or_eff,toplot,conditions,release);

kappa = k2/kreal-1;


