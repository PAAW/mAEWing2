function [KS,UB]=KS_opt_xAEWing3(g,gmax)

rho = 50; % From Martin's paper

ng=length(g);
% gmax = 0;


sum_constraints=0;

for num=1:ng
       
   sum_constraints=sum_constraints+exp(rho*(g(num)-gmax));
    
end

KS=1/rho*log(sum_constraints);

LB = gmax + log(ng)/rho - 2;
UB = gmax + log(ng)/rho;