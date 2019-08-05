% This function is to verify the SLP code for weight minimization of
% stiffened plate subjected to both buckling and stress constraints against
% NASTRAN optimization results.
% Contact: Wei Zhao 
% weizhao@vt.edu
%
%
clear all;
% close all;
clc;
addpath([pwd filesep 'Matlab-NASTRAN']);

addpath(genpath([ pwd filesep 'sensitivity_extract']));

addpath('sensitivity_extract\KS_function');
% addpath('fminslp');

addpath('slp_sqp');

%%

design_folder = [pwd filesep 'verification_stiffened_plate'];
main_folder = pwd;

nastran = 'C:\MSC.Software\MSC_Nastran\20170\bin\nast20170.exe ';

sol200_bdf_fname = 'dsoug4_only_sens_buckling';

DESVAR_labels = {'T-PLATE'
    'HATDIM2'};
%%
x0= [0.15
     1];

dxlb= [0.001
    0.001];
dxub=  [10.
    10.];

fit_function =@(x) compute_fitness_constraints_Canfield_SLP(x);
grad_function =@(x) compute_gradient_Canfield_SLP(x);


options.Display='iter';
options.MaxIter=100;
options.TolX =1e-6;
options.TolFun=1e-6;
options.TolCon=1e-6;
options.DerivativeCheck='on';
options.CompexStep='on';
% options.MoveLimit=1; % relative to X0; absolute, if options.TypicalX=1
%options.MoveReduction=0.8;
%options.TrustRegion='off';

%% SLP
options.TrustRegion='merit';
[x,f,exitflag,output]=slp_trust(fit_function,x0,options,dxlb,dxub,grad_function);
disp(' ')
disp('Final Design Variables, X')
disp(x)

% END
