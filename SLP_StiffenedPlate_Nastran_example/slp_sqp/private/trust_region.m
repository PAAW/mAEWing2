function TR=trust_region(TR,x,dx,vlb,vub,f0,f,gradf,g0,g,gradg,lambda)
% trust_region - determine trust ratio, move limits, binding constraints & penalties
%
%--Input/Output
%  TR.........   parameters and options for Trust Region (TR) strategy (structure)
%  .MoveLimit    Initial L-inf relative move limit for trust radius (default=0.2)
%  .MoveRelative Move limits relative to current x instead of x0
%  .Reject       Trust ratio value to reject new iterate    (default <0.1)
%  .Contract     Trust ratio value to contract trust region (default <0.2)
%  .Expand       Trust ratio value to expand   trust region (default >0.75)
%  .MoveExpand   Move Limit expansion factor              (default 2.0)
%  .MoveContract Move Limit contraction/reduction factor  (default 0.5)
%                for TrustRegion='off' apply every iter   (default 0.8)
%  .TrustRegion  Strategy options:
%                'simple' = exact L1 single penalty, mu*max(0,g) (default)
%                'merit'  = multiple-penalty parameter L1 merit function
%                'filter' = TR filter instead of merit descent function
%                'TRAM'   = TR adaptive move limit bounds (with merit)
%                'TRAM-filter'   combination
%                'simple-filter' combination
%
%--Input
%  x.......... design variable current values (vector)
%  dx......... change in design variable values from last iteration (vector)
%  vlb........ design variable lower bounds (vector)
%  vub........ design variable upper bounds (vector)
%  f0,f....... previous and current objective function values (scalar)
%  gradf...... gradient of objective function (vector)
%  g0,g....... previous and current constraint values (vector)
%  gradg...... gradients of constraints wrt design variables (matrix)
%  lambda..... Lagrange multipliers for constraints (structure)
%
%--Output
%
%  TR.delx....... bound on change in design variables (vector)
%  TR.rejected... flag for rejection or acceptance of current design iterate
%
%  Written by:    Robert A. Canfield
%  e-mail:        bob.canfield@vt.edu
%  Created:        10/ 2/16
%  Last modified:   4/24/18

%--Modifications
%  11/ 7/16 - Handle missing or null TypicalX
%  11/ 9/16 - Trust Region Approximation Method (TRAM)
%  11/11/16 - TR filter option
%  11/14/16 - simple exact L-1 penalty merit function
%  11/22/16 - TR.options
%   4/ 9/17 - Unblock filter & defaut is filter
%   4/18/17 - MoveRelative option added
%   4/23/17 - FCD modified for reducing constraint violation
%   5/29/17 - Filter block OK for feasible
%   1/24/18 - Adjust penalty/merit after accept/reject, not before
%   2/25/18 - smarter TypicalX
%   4/21/18 - account for TR increased approximate objective, numeric g=0
%             don't allow duplicate point through the filter
%   4/24/18 - TRAM line search only if FCD satisfied
%
% Wujek, B. A., and Renaud, J. E. "New Adaptive Move-Limit Management
% Strategy for Approximate Optimization, Part I," AIAA J. Vol. 36, No. 10,
% 1998, pp. 1911-1921
%
% Nocedal, J., and Wright, S.J. Numerical Optimization. New York: Springer, 2006.
% Algorithms 4.1 for TR & 15.1 for filter; Equation (15.4) for simple merit
%
% For the TR filter algorithm, also see 
% Fletcher, R., and Leyffer, S. 
% "Nonlinear programming without a penalty function," 
% Mathematical Programming Vol. 91, No. 2, 2002, pp. 239-269.
% doi: 10.1007/s101070100244

%% Process TR inputs and options
if isfield(TR,'TrustRegion'), TrustRegion=TR.TrustRegion; else, TrustRegion='on'; end
TR.trust = ~strcmpi(TrustRegion,'off');
TR.adapt = strncmpi(TrustRegion,'TRAM',4);
TR.SimpleTrust = contains(lower(TrustRegion),'simple');
Filter = (strcmpi(TrustRegion,'on') ||  contains(lower(TrustRegion),'filter'))...
                                    && ~contains(lower(TrustRegion),'merit');
if TR.adapt && isempty(which('trust_adapt','in','trust_region'))
   error('trust_region:noTRAM','Trust Region Adaptive Move Limit not available, yet.')
end
if ~isfield(TR,'Reject')    || isempty(TR.Reject),    TR.Reject    = 0.10; end
if ~isfield(TR,'MoveLimit') || isempty(TR.MoveLimit), TR.MoveLimit = 0.20; end
if ~isfield(TR,'Contract')  || isempty(TR.Contract),  TR.Contract  = 0.20; end
if ~isfield(TR,'Expand')    || isempty(TR.Expand),    TR.Expand    = 0.75; end
if ~isfield(TR,'MoveExpand')|| isempty(TR.MoveExpand),TR.MoveExpand= 2.0;  end
if ~isfield(TR,'MoveRelative') || isempty(TR.MoveRelative),TR.MoveRelative='off'; end
if ~isfield(TR,'MoveContract') || isempty(TR.MoveContract)
   if isfield(TR,'MoveReduce') % backward compatibility
      TR.MoveContract=TR.MoveReduce;
   elseif isfield(TR,'MoveReduction')
      TR.MoveContract=TR.MoveReduction;
   elseif TR.trust
      TR.MoveContract = 0.5;
   else
      TR.MoveContract = 0.9;
   end
end
MoveRelative = strcmpi(TR.MoveRelative,'on');
% Other TR fields used:
%       TR.cutoff
%       TR.SimpleTrust
%       TR.options.TolCon
%       TR.options.TolX
%       TR.options.TypicalX

%% Local variables
[mg0,mj0] = max(g0);
[mg, mj]  = max(g);

%% Initialization (before first optimization iteration)
if ~isfield(TR,'delx') || isempty(TR.delx)
   if ~isfield(TR,'TolFun'), TR.TolFun=TR.options.TolFun; end
   if ~isfield(TR,'TolCon') || isempty(TR.TolCon)
                             TR.TolCon=TR.options.TolCon; end
   if ~isfield(TR,'TolX'),   TR.TolX  =TR.options.TolX;   end
   if isfield(TR,'TypicalX') && ~isempty(TR.TypicalX) && ~ischar(TR.TypicalX)
      if length(TR.TypicalX)~=length(x)
         error('trust_region:TypicalX','TR.TypicalX not length of x')
      else
         TypicalX=TR.TypicalX(:);
      end
   else
      if all( abs(x) > TR.TolX )
         TypicalX = x(:);
      else
         CenterX = abs(vub + vlb)/2;
         if all( CenterX > TR.TolX ) && all( CenterX < 1e8 )
            TypicalX = CenterX(:);
         else
            TypicalX = ones(size(x(:)));
         end
      end
   end
   TR.delx = max(TR.options.TolX, TR.MoveLimit*abs(TypicalX));
   if isinf(g)
      Lambda = 0;
   else
      Lambda = norm(gradf) ./ max(eps,sqrt(diag(gradg'*gradg)));
   end
   if TR.SimpleTrust
      TR.penalty = max(Lambda);
   else
      TR.penalty = Lambda;
   end
   TR.Merit = trust_merit(f,g,TR.penalty);
   TR.bindLM = find(g>=min(0,mg)-TR.cutoff);
   TR.contracted = false;
   TR.rejected=[];
   TR.output.RadiusFraction=[];
   TR.filter.flag = Filter;
   TR.filter.f  = f0;
   TR.filter.mg = mg0;
   TR.filter.mj = mj0;
   return
else
%% Determine binding constraints and their penalties
   tol = optimget(optimset('quadprog'),'TolCon'); % tolerance for active constraints
   bound = ( (lambda.lower>tol & x>vlb+tol) ...
           | (lambda.upper>tol & x<vub-tol) ) & (TR.delx-abs(dx))<tol;
   % Recover saved values from TR
   MoveLimit0 = TR.MoveLimit;
   Merit0     = TR.Merit;
   penalty0   = TR.penalty;
   delx0      = TR.delx;
end

%% Determine binding constraints and their penalties
TR.Bound = any( bound ); % move-limit-bound not at variable bounds
if TR.stat<0 % linprog did not converge; Lagrange multipliers unreliable
   bindLM=[];
else
   bindLM = find(g>=min(0,mg)-TR.TolCon ...
      | lambda.ineqlin>tol*max(lambda.ineqlin));
   TR.bindLM = bindLM;
end

%% Actual and approximate merit function values
f1 = TR.fapx;
g1 = TR.gapx;
[Merit,penalty] = trust_merit( f, g, penalty0, lambda, bindLM );
MeritApx        = trust_merit( f1,g1,penalty );
Violated     = mg > TR.TolCon;
MoreViolated = mg > mg0 && Violated;

%% Objective/Constraint-Violation Filter
if all( f < TR.filter.f & max(0,mg) <= max(0,TR.filter.mg) )
   FilterStr = '*'; % Dominates all previous designs
   NonDom = 2;
elseif any(TR.filter.f <= f+eps & max(0,TR.filter.mg) <= max(0,mg+eps))
   FilterStr = '-'; % Dominated
   NonDom = 0;
else
   FilterStr = '+'; % Nondominated
   NonDom = 1;
end

%% Trust Region Strategy
if TR.trust
   % Trust Ratio
   tolg = TR.TolCon/10;
   j = unique([mj0,mj,find(g0>tolg | g>tolg)']);
   g0j = max(g0(j),tolg);
   g1j = max(g1(j),tolg);
   if f1 < f0 % approximate objective decreased (satisfies FCD wrt eps)
      TrustRatio_f = (f0-f) / max(f0-f1,eps);
   else       % approximate objective increased (maybe to make feasible)
      den = max(abs(f0-f),eps)*sign(f0-f);
      TrustRatio_f = (f0-f1) / den;
   end
   TrustRatio_g = min( (g0(j)-g(j)) ./ (g0(j)-g1(j)) );
   TrustRatio_M = (Merit0-Merit) / (Merit0-MeritApx);
   if isfield(TR,'MPEA') && strcmpi(TR.MPEA,'beta')
      TR.Ratio_f = TrustRatio_f;
      TR.Ratio_g = zeros(size(g));
      TR.Ratio_g(j) = (g0j - g(j)) ./ max(abs(g0j - g1j),eps);
   end
   if Filter
      if all([g0j,g1j] <= tolg)
         [TR.Ratio,r] = min([TrustRatio_f, inf, TrustRatio_M]);
      else
         [TR.Ratio,r] = min([TrustRatio_f, TrustRatio_g]);
      end
      fgmstr='fgm';
      FilterStr=[fgmstr(r),' ',FilterStr];
   else
      TR.Ratio = TrustRatio_M;
   end

%% Accept or reject new design point
   % Increase penalty if violation increased
   if Merit<Merit0 && MoreViolated && ~TR.SimpleTrust
     penalty(mj) = penalty(mj) + (Merit0-Merit)/mg;
     Merit = Merit0; % Reset merit to force contraction
   end
   if Filter
      TR.rejected = ~NonDom && (MoreViolated || ~Violated); % Dominated (unblock ~MoreViolated)
   elseif TR.SimpleTrust
      TR.rejected = TR.Ratio < TR.Reject; % trust ratio too low
   else
      TR.rejected = (f > f0 && max(0,mg) >= max(0,mg0)) || ... % dominated
                    (Merit > Merit0 && MoreViolated); % merit & mg too high
   end

   mgapx = max(TR.gapx);
   %%    Fraction of Cauchy Decrease
   FCD = MeritApx < Merit0 ... % Sufficient FCD
         || (mgapx > TR.TolCon && mgapx < mg0);
   if TR.rejected
      FCDstr=' Rejected';
      Merit = Merit0;
   else
      FCDstr=[];
   end
   if ~FCD, FCDstr=[FCDstr ' ~FCD']; end

   %% Filter or Merit function determines Trust Region
   Expand = TR.Ratio > TR.Expand && ~MoreViolated && ~TR.contracted && TR.Bound;
   if Filter
      if f>f0 && max(0,mg)>max(0,mg0) % Dominated by last point
         NonDomLast = 0;
      elseif f<f0-TR.TolFun && mg<max(TR.TolCon,mg0-TR.TolCon) % Dominates last point
         NonDomLast = 2;
      else
         NonDomLast = 1 ;
      end
      Contract = TR.rejected || (NonDomLast<2 && TR.Ratio<TR.Contract);
      Expand   = NonDom >= 1 && Expand;
   else
      Contract = TR.rejected || ~FCD || TR.Ratio < TR.Contract;
   end

   %% Adapt Move Limit bounds (Trust Region Approximate Model)
   MoveReduce = TR.MoveContract;
   MoveExpand = TR.MoveExpand;
   TRAMstr=[];
   if TR.adapt && FCD
      TR.Merit      = Merit;
      TR.Merit0     = Merit0;
      [R,TRAMstr] = trust_adapt(f0,f,f1,g0,g,g1,gradf,gradg,dx,penalty0,penalty,TR);
      MoveReduce = R;
      MoveExpand = R;
   end
   if MoveReduce >= 1 && Contract
      error('slp_trust:MovRed','Trust Region Strategy Move Reduction must be <1')
   end

   %% Apply Move Limit factor to design variable change (delx)
   if Contract
      TR.delx      = MoveReduce*TR.delx;
      TR.radius    = MoveReduce;
      if all(TR.delx < TR.TolX)
         warning('trust_region:TolX','Move limits contract < TolX')
      end
   elseif Expand && TR.Bound
      TR.delx(bound) = MoveExpand*TR.delx(bound);
      TR.radius      = MoveExpand;
   else
      TR.radius = 1;
   end
   TR.MoveLimit  = TR.radius*MoveLimit0;
   TR.contracted = TR.MoveLimit < MoveLimit0;
   TR.filter.f(end+1)  = f;
   TR.filter.mg(end+1) = mg;
   TR.filter.mj(end+1) = mj;
   TR.filter.ftype = Merit < Merit0 && FCD;
   TR.filter.str = strcat(FilterStr,TRAMstr,FCDstr);

else
   %% No trust region strategy (ad hoc move limits)
   if TR.Bound || mg>=mg0
      if MoveRelative
         %     TR.delx = max(TR.TolX, TR.MoveLimit*abs(x(:)));
         TR.delx = TR.MoveLimit*abs(x(:));
      else
         TR.delx = TR.MoveContract*TR.delx;
      end
   end
   TR.MoveLimit = TR.MoveContract*MoveLimit0;
   TR.radius = [];
   TR.Ratio  = [];
   TR.rejected=false;
   TR.filter.str = FilterStr;
end

%% Save local variables to return for next iteration
TR.Merit      = Merit;
TR.Merit0     = Merit0;
TR.bound      = bound;
TR.penalty    = penalty;
if ~isempty(TR.radius)
   TR.output.RadiusFraction(end+1) = min(TR.delx./delx0)/TR.radius;
end

end