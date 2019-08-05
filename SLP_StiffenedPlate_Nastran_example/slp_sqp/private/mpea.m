function [y,dylb,dyub,gradfy,gradgy,TR]=mpea(TR,x,dxlb,dxub,f,g,gradf,gradg,...
                                        lambda0,x0,xp,xr,g0,fp,gp,gradfxp,gradgxp)
%  Multipoint Exponential Approximation Transformation
%
% Canfield, R. A. "Quadratic Multipoint Exponential Approximation:
% Surrogate Model for Large-Scale Optimization," 12th World Congress of
% Structural and Multidisciplinary Optimization (WCSMO12). Vol. Advances in
% Structural and Multidisciplinary Optimization, Springer International
% Publishing, Braunschweig, Germany, 2017, pp. 648-661.
%
%  Written by:    Robert A. Canfield
%  e-mail:        bob.canfield@vt.edu
%
%  Created:        4/2/17
%  Last modified:  4/7/19

%--Innpt
%  lambda0... Lagrange multpliers from last iteration
%  x0........ most recent point (may be rejected)
%  xp........ previous point with gradients
%  xr........ rejected point (empty if accepted)
%  g0........ most recent constraint values (may be rejected w/o gradients)

%--Modifications
%   2/1/18 - variable trust ratios govern exponent instead of statistics
%   2/2/18 - wash out memory with historical trust ratios
%   2/8/18 - refactored as mpea and qmea
%  4/19/18 - option to limit to current sensitivity
%  4/22/18 - check function matching
%  4/24/18 - dylb/ub check for crossing zero
%   5/3/18 - save internal variables instead of appending to TR
%   4/7/19 - reset Phis after slack variable from last iter

% Previous Lagrangian sensitivity for legacy MPEA
persistent g_linear iter_mpea Phis Pmean Psigma pmask

%% Return inverse transformation 
if nargin==2 && nargout==1 % really x=y^(-p)
   if isfield(TR,'p') && ~isempty(TR.p) && ~all(TR.p==1)
      y = x.^(1./TR.p);
   else
      y = x;
   end
   return
end

%% No previous iteration or linear, return y=x
if isempty(gradfxp)
   y    = x;
   dylb = dxlb;
   dyub = dxub;
   gradfy = gradf;
   gradgy = gradg;
   g_linear=true(size(g));
   iter_mpea = 0;
   TR.p=ones(size(x));
else % there was previous point and not linear (all p=1)

%% Lagrange Multipliers before and after current iteration
if isfield(TR,'Lambda')
   LM = max( lambda0.ineqlin, TR.Lambda );
else
   LM = lambda0.ineqlin;
end

%% Identify linear constraints and p=1 exponents
[ndv,ncon] = size(gradg);
% [mg,mj] = max(g);
if isempty(xr)
   dx = x0 - x; % Expansion from current to previous accpeted point
else
   dx = xr - x; % Expansion from current to previous rejected point
   x0 = xr;
end
gapx = g + gradg.'*dx;
g_linear = abs(g0-gapx) < TR.TolCon & g_linear; % Linear constraints
if ~isempty(gradgxp)
   g_linear = g_linear & all( abs(gradg-gradgxp) < TR.TolCon, 1).';
end
a_lin = g_linear & LM > TR.TolCon & max(g,g0) > -TR.TolCon; % active linear g
i_bnd = TR.bound & abs(dx)<eps ...
   & (lambda0.lower | lambda0.upper);  % x_i at bound
i_fix = i_bnd | any(gradg(:,a_lin),2) ... % active linear constraint
   | (x+dxlb) < 0 ...           % sign change of x not allowed
   | abs(dx) < eps | x.*x0 < 0; % p=1 if unchanged variables or sign change

%% Update exponents if last point accepted
if TR.rejected % No new gradient to update p
   p = TR.p;
   p(i_fix) = 1;
else % accepted point
   p = ones(size(x));
   if isempty(Phis) || length(Phis)~=length(gradf)
      Phis = zeros(size(gradf));
      TR.Pj = zeros(size(gradf));
      Pmean = zeros(size(gradf));
      Psigma = zeros(size(gradf));
   else
      for i=any(gradg(:,a_lin),2)
         P_lin = gradg(i,a_lin).*LM(a_lin)';
         Phis(i,a_lin) = P_lin; % store sensitivity
         [~,TR.pj(i)] = max(abs(P_lin)); % record linear constraint index
      end
   end
   
   %% Current constraint Lagrangian sensitivity only
   if strcmpi(TR.MPEA,'current') % simple approach
      % Constraint gradient sensitivity
      P = repmat(LM',numel(x),1).*gradg;
      lset = abs(P) > TR.TolFun & gradg.*gradgxp>0 ...
           & repmat(gradf,1,size(gradg,2)).*gradg < 0;
      dx = xp - x;
      if isfield(TR,'qmea') && isfield(TR.qmea,'Beta') && any(TR.qmea.Beta)
         G_beta_d = TR.qmea.basis*diag(TR.qmea.Beta)*TR.qmea.basis'*dx;
         G_Beta_D = repmat(G_beta_d,1,size(gradg,2));
      else
         G_beta_d = false(size(gradf));
         G_Beta_D = zeros(size(gradg));
      end
      XbyX = repmat(log(xp./x),1,ncon);
      pij = zeros(size(gradg));
      pij(lset) = 1 + log(max(eps,(gradgxp(lset)-G_Beta_D(lset)) ./ gradg(lset)))...
                    ./XbyX(lset);
      for i=find(~i_fix & any(lset,2))'
         if ~any( P(i,a_lin) )      % Leave critical linear constraints p=1
            [~,j]=max(abs(P(i,:))); % Dominant nonlinear constraint
            TR.pj(i) = j;           % record constraint index
            p(i) = pij(i,j);
            % Check stability of recursive exponent eq.
            if (G_beta_d(i))
               phi = @(Gbd) 1 + log(max(eps,(gradgxp(i,j)-Gbd) ./ gradg(i,j)))...
                  ./XbyX(i,j);
               dphi = imag( phi( G_Beta_D(i,j) + 1i*eps )) / eps;
               if (phi(G_Beta_D(i,j)) ~= p(i))
                  error('mpea:p','p(i) is not same as phi(i)')
               elseif abs(dphi) >= 1
                  p(i) = phi(0);
               end
            end
         end
      end
      i = ~any(lset,2) & ~i_fix ... % rest from gradf
        & gradf.*gradfxp > eps; % consistent sign of objective derivative
      if any(i)
         phi = @(Gbd) 1 + log(max(eps,(gradfxp(i)-Gbd) ./ gradf(i)))...
            ./ (log(x0(i)./x(i)));
         p(i) = phi( G_beta_d(i) );
         % Check stability of recursive exponent eq.
         if (G_beta_d(i))
            phi = @(Gbd) 1 + log(max(eps,(gradfxp(i)-Gbd) ./ gradf(i)))...
                ./ (log(x0(i)./x(i)));
            dphi(i) = imag( phi( G_beta_d(i) + 1i*eps )) / eps;
            for j = i & abs(dphi(i)) >= 1
               p(j) = phi(0);
            end
         end

      end
      
   elseif strcmpi(TR.MPEA,'beta') % beta-development trial version
      %% Weight Lagrangian sensitivities by trust ratios
      %  Recreate MPEA using current gradients and previous exponents
      [gradfy,gradgy,y]=transf_x2y(TR.p,gradf,gradg,x);
      dy = x0.^TR.p - y;
      fapx0 = f + gradfy.'*dy;
      gapx0 = g + gradgy.'*dy;
   
      % Wash out prior Lagrangian sensitivity with latest trust ratios
      i = TR.Pj==0;
      TR.P(i) = max(0,min(1,TR.Ratio_f)) * TR.P(i);
      i = TR.Pj>0;
      j = TR.Pj(i);
      TR.P(i)  = max(0,min(1,TR.Ratio_g(j))) .* TR.P(i);

      % Trust ratio for each intermediate variable in objective approximation
      fapxi = fapx0 + gradfy(i).*((x0(i).^p(i)-x(i).^p(i)) - dy(i));
      tr_f = zeros(size(gradf));
      tr_f(i) = max(0, min(1, (fp - f) ./ (fapxi - f)));
      P = abs(tr_f.*gradf); % objective sensitivity weighted by trust ratio
      i = i & P > Phis;
      Phis(i) = P(i); % replace Lagrangian sensitivity with larger weighted obj sens.
      TR.p(i) = p(i);
      TR.Pj(i) = 0;
      i_lin = any(gradg(:,a_lin),2);
      [~,TR.Pj(i_lin)]=max( abs(gradg(i_lin,:)), [], 2 );
      
      % Dominant constraint derivatives determine exponent
      lset = gradg.*gradgxp>0; % & repmat(gradf,1,size(gradg,2)).*gradg < 0;
      lset(i_fix,:) = false;
      lset(:,g < -TR.TolCon) = false;
      XbyX = repmat(log(x0./x),1,ncon);
      pij = zeros(size(gradg));
      pij(lset) = 1 + log(gradgxp(lset) ./ gradg(lset))./XbyX(lset);
      
      % Trust ratio for each intermediate variable in constraint approximation
      DY0 = repmat(dy,1,ncon);
      DY  = repmat(x0,1,ncon).^pij - repmat(x, 1,ncon).^pij;
      gapxi = repmat(gapx0.',ndv,1) + gradgy.*(DY - DY0);
      tr_g = max(0, min(1, repmat((gp - g)',ndv,1) ./ (gapxi - repmat(g',ndv,1))));
      LMgradg = abs(repmat(LM',ndv,1).*gradg);
      P = zeros(size(gradg));
      P(lset) = abs(tr_g(lset).*LMgradg(lset));
      [P,j] = max(P,[],2); % collapse matrix of Lagrangian sensitivity to one per var.
      i = find( P > Phis );
      Phis(i) = P(i); % replace Lagrangian sensitivity with larger weighted sens.
      TR.p(i) = pij( sub2ind(size(pij), i, j(i)) ); % update associated exponents
      TR.Pj(i) = j(i);
      p = TR.p;
      clear dy DY DY0 fapxi gapxi LMgradg
      
   else % legacy approach 'on'
      %%  Scan Lagrangian sensitivity table to determine exponent
      % Objective gradient governs exponent if constraint gradients won't
      i = ~i_fix & gradf.*gradfxp  > eps; % consistent sign of objective derivative
      p(i) = 1 + log(gradfxp(i)./gradf(i)) ./ (log(x0(i)./x(i)));
      P = abs(gradf);
      i = i & P > Phis;
      Phis(i) = P(i); % replace Lagrangian sensitivity with larger weighted obj sens.
      TR.p(i) = p(i); % store variable exponent
      TR.Pj(i) = 0;   % record that objective (0 index) governed choice of p
      p(~i) = TR.p(~i); % recover previous exponents for other variables
      
      % Constraint gradient sensitivity
      P = repmat(LM',numel(x),1).*gradg;
      Pmean  = mean(P,2);
      lastP = iter_mpea>1;
      if lastP
         Pmean  = ((iter_mpea-1)*Pmean +Pmean)/iter_mpea;
         Psigma = std([(iter_mpea-1)*P, Pmean]/iter_mpea,0,2);
      else
         Psigma = std(P,0,2);
         pmask = false(size(p));
      end
      Pupper = min(P, Pmean + Psigma);
      Plower = max(P, Pmean - Psigma);
      lset = (P >= Pupper | P <= Plower) & abs(P) > TR.TolFun...
         & gradg.*gradgxp>0 & repmat(gradf,1,size(gradg,2)).*gradg < 0;
      
      % All constraint gradients in Lagrangian used if no dominant constraints
      i = ~i_fix & gradg.*gradgxp*LM > eps & ~any(lset,2);
      p(i) = 1 + log(gradg(i,:)*LM./(gradgxp(i,:)*LM))...
         ./ (log(x(i)./xp(i)));
      
      % Critical linear constraints
      p(any(abs(P(:,a_lin)) > TR.TolFun,2)) = 1;
      i_lin = any(gradg(:,a_lin),2);
      [~,TR.Pj(i_lin)]=max( abs(gradg(i_lin,:)), [], 2 );
      
      % Dominant nonlinear constraint
      XbyX = repmat(log(xp./x),1,ncon);
      pij = zeros(size(gradg));
      pij(lset) = 1 + log(gradgxp(lset) ./ gradg(lset))./XbyX(lset);
      for i=find(~i_fix & any(lset,2))'
         if ~any( P(i,a_lin) ) % Leave critical linear constraints p=1
            if pmask(i)     % include previous p, if set by dominant P
               pijk = [pij(i,lset(i,:)), TR.p(i)];
               lsetj = [find(lset(i,:)), TR.Pj(i)];
            else
               pijk = pij(i,lset(i,:));
               lsetj = find(lset(i,:));
            end
            [pm1,km1]=min(abs(pijk-1));
            [pp1,kp1]=min(abs(pijk+1));
            if pp1 < pm1
               p(i) = pijk(kp1);
               TR.Pj(i) = lsetj(kp1); % record governing constraint index
            else
               p(i) = pijk(km1);
               TR.Pj(i) = lsetj(km1); % record governing constraint index
            end
            % 	      [~,j]=max(abs(P(i,:)));
            % 	      p(i) = pij(i,j);
            pmask(i)=true; % an exponent was chosen by constraint
         end
      end
   end
   clear XbyX
   if any(isnan(p(:)) | isinf(p(:)))
      disp('mpea: p=NaN or inf reset to p=+/-4')
   elseif any(p(i_fix)~=1)
      warning('mpea:i_fix','Fixed variables has p~=1')
   end
   p=min(4,max(-4,p));
end % accepted point

%% MPEA variable transformation
[gradfy,gradgy,y,dylb,dyub]=transf_x2y(p,gradf,gradg,x,dxlb,dxub,TR.TolX);
TR.p = p;
if TR.rejected
   return
end

end % of there was (not) a previous point

%% Save current iterate data for trust ratio and next iteration
%  if not iterating on beta
if ~isfield(TR,'QMEA') || ~isfield(TR,'qmea') || ~isfield(TR.qmea,'Beta')
   % Save latest grad for p-beta iteration
   if isfield(TR,'QMEA') && ~strcmpi(TR.QMEA,'off')
      TR.gradfyp=gradfy;
      TR.gradgyp=gradgy;
   end
end
end % mpea function



%% Transform x to y sub-function
   function [gradfy,gradgy,y,dylb,dyub]=transf_x2y(p,gradf,gradg,x,dxlb,dxub,TolX)
      %% Transform to intermediate variables
      if all( p==1 )
         y = x;
         if nargout>3
            dylb = dxlb;
            dyub = dxub;
         end
         gradfy = gradf;
         gradgy = gradg;
      else
         dxdy = x.^(1-p)./p;
         gradfy = gradf.*dxdy;
         gradgy = repmat(dxdy,1,size(gradg,2)).*gradg;
         y = x.^p;
         if nargout>3
            dylb = dxlb;
            dyub = dxub;
            i = p>0 & p~=1;
            dylb(i) = (max(TolX,x(i)+dxlb(i))).^p(i) - y(i);
            dyub(i) =          (x(i)+dxub(i) ).^p(i) - y(i);
            i = p<0;
            dylb(i) =           (x(i)+dxub(i)).^p(i) - y(i);
            dyub(i) = (max(TolX,x(i)+dxlb(i))).^p(i) - y(i);
            dyub    = max( dyub, TolX.^min(1,p) );
            if any(dylb>dyub),     error('mpea:dylub','wrong bound'), end
            if any(~isreal(dylb)), error('mpea:dylb','complex dylb'), end
            if any(~isreal(dyub)), error('mpea:dyub','complex dyub'), end
         end
      end
   end