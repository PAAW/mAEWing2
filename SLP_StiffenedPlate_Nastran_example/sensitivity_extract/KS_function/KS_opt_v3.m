function [KS,UB]=KS_opt_v3(g,gmax,rho)

% K-S constraints

% gi/maxvalue<=1



% for too many constraints, use first 1000 instead

% if length(g)>=1000
%     
%    for constraint_num=1:length(g)
%        
%       diff_2_max(constraint_num)=g(constraint_num)-maxvalue;
%        
%        
%    end
%     
%    
%    
%     
%     
% end


% rho=50;length(responses);

ng=length(g);


sum_constraints=0;

for num=1:ng
       
   sum_constraints=sum_constraints+exp(rho*(g(num)-gmax));
    
end

KS=1/rho*log(sum_constraints);

LB = gmax + log(ng)/rho - 2;
UB = gmax + log(ng)/rho;