%% runBeamGV
% Script to Run Gary Vanderplaats cantilever Beam with SQP and SLP. 
% N beam segments and 2N design variables, taken from 
% Vanderplaats (1984) Example 5-1, pp. 147-150.

%% Initialize variables
clear; clc
options.Display = 'iter';
options.MaxIter = 40;
options.TolX   = 0.5;
options.TolFun = 0.5;
options.TolCon = 1e-4;
options.MoveLimit = 0.2;

%% Small and Large-Scale number of variables
for N = [5 200] % number of beam segments
   x0 = [5*ones(1,N), 40*ones(1,N)];
   xlb = [ ones(1,N),  5*ones(1,N)];
   options.TypicalX = x0;
   disp(' ')
   disp(['N=',num2str(N),' beam segments'])

%% SQP with BFGS
   disp('SQP with BFGS'), tic
   [dv2,out2,lam2] = sqp(@fbeamGV,x0,options,xlb,[],@gbeamGV); toc
   disp(out2)

%% SQP with Exact Hessian
   disp('SQP with exact Hessians')
   options.HessFun=@HbeamGV; tic
   [dv3,out3,lam3] = sqp(@fbeamGV,x0,options,xlb,[],@gbeamGV); toc
   disp(out3)
   options=rmfield(options,'HessFun');

%% SLP Trust
%  with and without active set strategy
   if N>50
      options.TrustRegion='TRAMfilter'; 
      tic
      [~,fa,stat,outa] = slp_trust(@fbeamGVa,x0,options,xlb,[],@gbeamGVa), toc
   else
      options.TrustRegion='filter'; tic
      [~,f1,sta1,out1] = slp_trust(@fbeamGV,x0,options,xlb,[],@gbeamGV)
   end

%% fmincon with analytic gradients
   A = [diag(-20.*ones(N,1)),diag(ones(N,1))]; % Linear Constraints
   b_U     = zeros(N,1);                       % Bound on linear constraints
   disp(' ')
   disp('fmincon fails with looser tolerances used for SQP and SLP_Trust')
   Options=optimoptions(@fmincon,'Display','iter','GradObj','on','GradConstr','on',...
                                 'TolX',0.5,'TolFun',0.5,'TolCon',1e-4);
   [~,~,flag,out]=fmincon(@GVbeam_obj,x0,A,b_U,[],[],xlb,[],...
                          @GVbeam_nlcon,Options)
   Options=optimoptions(Options,'Algorithm','sqp');
   [~,~,flag,out]=fmincon(@GVbeam_obj,x0,A,b_U,[],[],xlb,[],...
                          @GVbeam_nlcon,Options)

 %% fmincon, tight tolerances
    disp('fmincon with tighter default tolerances')
    Options=optimoptions(@fmincon,'Display','iter','GradObj','on','GradConstr','on');
    tic
    [~,~,flag,out]=fmincon(@GVbeam_obj,x0,A,b_U,[],[],xlb,[],...
       @GVbeam_nlcon,Options)%#ok<*ASGLU,*NOPTS>
    toc
%     disp('fmincon treating all constraints as nonlinear')
%     tic
%     [~,~,flag,out]=fmincon(@GVbeam_obj,x0,[],[],[],[],xlb,[],...
%        @GVbeam_con,Options);
%     toc
    Options=optimoptions(Options,'Algorithm','sqp'); tic
    [~,~,flag,out]=fmincon(@GVbeam_obj,x0,A,b_U,[],[],xlb,[],...
                           @GVbeam_nlcon,Options), toc

end