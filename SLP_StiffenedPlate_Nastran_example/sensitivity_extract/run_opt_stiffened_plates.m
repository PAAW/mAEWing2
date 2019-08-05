
clear all;
clc;
addpath('C:\Users\weizhao\Documents\Matlab-NASTRAN');

filename = 'C:\Users\weizhao\Documents\PAAW\N+3\MATLAB_NASTRAN_Opt\verification_stiffened_plate\dsoug4_only_sens.f06';


addpath([pwd filesep 'KS_function\'])


%%

options = optimoptions('fmincon','SpecifyObjectiveGradient',true);
options.SpecifyConstraintGradient=true;
% options.ConstraintTolerance = 1e-3;

options.MaxFunctionEvaluations = 1e4;



x0 =[  1.1271E-01 8.3742E-01 
];
lb=[0.001    0.001];

ub = [10.    10.];


[x,fval,exitflag,output] = fmincon(@fitness,x0,[],[],[],[],lb,ub,@nonlinear_constraints_SP_just_disp,options);
%% Displacement