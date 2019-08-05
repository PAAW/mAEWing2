%% run Simple explict example

%% Initialize variables
clear; clc
options.Display = 'iter';
options.MaxIter = 50;
options.TolX   = 0.001;
options.TolFun = 0.1;
options.TolCon = 1e-4;
options.OptimalityTolerance = 0.1;
options.MoveLimit = 0.2;
options.ComplexStep = 'on';

x0 = ones(1,3);
xlb = 0.1*x0;
xub =  10*x0;

%% SLP Trust
options.TrustRegion='filter'; options.MPEA='off';
[x1,f1,stat1,outslp,LMslp] = slp_trust(@fsimple,x0,options,xlb,xub)

%% MPEA legacy
options.TrustRegion='filter'; options.MPEA='legacy'; options.QMEA='off';
[x2,f2,stat2,outmpea] = slp_trust(@fsimple,x0,options,xlb,xub)

%% MPEA Trust
options.TrustRegion='filter'; options.MPEA='on'; options.QMEA='off';
[x3,f3,stat3,outmpea0] = slp_trust(@fsimple,x0,options,xlb,xub)

%% QMEA Trust (quadprog)
options.MPEA='on'; options.QMEA='on';
[x4,f4,stat4,outqmea] = sqp_trust(@fsimple,x0,options,xlb,xub)