    function M = KindOfDynamics( NameOfModel )
% KINDOFDYNAMICS Summary of this function goes here
% Detailed explanation goes here
% % Models of dynamics
global A a h
switch NameOfModel
    case 'MAK'
        % Biochemical
        M0 = @(xi) 1-xi;
        M1 = @(xi) -xi;
        M2 = @(xj) xj;
        M = { M0, M1, M2};
    case 'PD'
        % Birth-Death Processes
        M0 = @(xi) -xi.^2;
        M1 = @(xi) 1;
        M2 = @(xj) xj.^3; 
        M = { M0, M1, M2};        
    case 'MM'
        % Gene Regulatory
%         a=1; h=2;
        M0 = @(xi) -xi.^a;
        M1 = @(xi) 1;
        M2 = @(xj) xj.^h./(1+xj.^h);
        Ginv = @(x) ( (x./(1-min(x,1))).^(1/h) ) ;
        Rinv = @(x) x.^(1/a);
        M = { M0, M1, M2 ,Ginv, Rinv};
    case 'SIS'
        % epidemics
        beta = 1;
        M0 = @(xi) -beta*xi;
        M1 = @(xi) 1-xi;
        M2 = @(xj) xj;
        M = { M0, M1, M2};        
    case 'Eco'
        % Ecology
        % Allee + migration        
        M0 = @(xi) 2 -   xi.*(xi-3).*(xi-6);
        M1 = @(xi) xi;
        M2 = @(xj) xj; %./(1+xj);
        Ginv = @(x) x;
        M = { M0, M1, M2, Ginv};        
    case 'Voter'
        % Voter model
    case 'opinions'
        M0 = @(xi) -(xi);
        M1 = @(xi) 1;
        M2 = @(xj) tanh(xj);
        Ginv = @(x) atanh(x);
        M = { M0, M1, M2, Ginv}; 
    case 'Simple'
        % simple invented example
        M0 = @(xi) -xi.^3;
        M1 = @(xi) 1;
        M2 = @(xj) xj;
        Ginv = @(x) x;
        M = { M0, M1, M2, Ginv};                
    case 'Glauber'
        M0 = @(xi) -(xi);
        M1 = @(xi) 1;
        M2 = @(xj) tanh(xj);
        Ginv = @(x) atanh(x);
        M = { M0, M1, M2, Ginv};                
    case 'Ising-Sch'      
        M0 = @(xi) -xi-xi.^3;
        M1 = @(xi) 1;
        M2 = @(xj) xj;
        Ginv = @(x) x;
        M = { M0, M1, M2, Ginv};                
    case 'Neural'
        mu = 5;
        M0 = @(xi) -xi;
        M1 = @(xi) 1;
        M2 = @(xj) 1./(1+exp(-xj+mu));
        Ginv = @(x) atanh(2*x-1)+1/2;
        M = { M0, M1, M2, Ginv};              
    otherwise
        disp('unknown model')
end
