function rhoc = find_rho_c(n,M,conditions,release,holding_value,x_th)
global x0
a = 0;
b = n;
while b-a>1
   m = floor(mean([a,b]));
   x0 = ( (1:n)<=m )' * holding_value;   
   xm = SolveOdes(x0,M,'eff',0,conditions,release);
   if xm<x_th
       a = m;
   else
       b = m;
   end
end

rhoc = b/n;