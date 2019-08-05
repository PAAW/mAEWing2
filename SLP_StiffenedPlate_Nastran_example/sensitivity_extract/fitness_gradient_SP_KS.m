function [f,gradf,c,gradc] = fitness_gradient_SP_KS(x0,...
    design_folder,main_folder,nastran,sol200_bdf_fname,DESVAR_labels)


% save x0-fitness.txt x0 -ascii -append

%% WRITE DESVAR.BDF


updated_DESVAR(x0);


%% RUN NASTRAN SOL 200

% % design_folder = 'C:\Users\weizhao\Documents\PAAW\N+3\MATLAB_NASTRAN_Opt\verification_stiffened_plate';
% % main_folder = 'C:\Users\weizhao\Documents\PAAW\N+3\MATLAB_NASTRAN_Opt\sensitivity_extract';
% % 
% % nastran = 'C:\MSC.Software\MSC_Nastran\20170\bin\nast20170.exe ';
% % 

% % sol200_bdf_fname='dsoug4_only_sens';
cd(design_folder);
call_NASTRAN_sol200=[nastran ' ' design_folder filesep sol200_bdf_fname '.bdf' '  memory=2g'];
system(call_NASTRAN_sol200);
!del *.xdb
!del *.1
% check program end or not
if isunix
    fname = [design_folder filesep sol200_bdf_fname];
    program_nastran_determine(fname)
end
%
cd(main_folder);


%%
% % DESVAR_labels = {'T-PLATE'
% %     'HATDIM2'};
filename  = [design_folder filesep sol200_bdf_fname '.f06'];

[initial_values,senstivity_matrix] = read_general_sensitivity(filename,DESVAR_labels,15,'WEIGHT');

f = initial_values;

gradf = senstivity_matrix;


% save weight_iteration.txt f -ascii -append

%% gradient

%% disp

C_LB = 0.1; %should provide positive value
C_UB = 0.1;

[initial_values,senstivity_matrix] = read_general_sensitivity(filename,DESVAR_labels,13,'DISP');

c(1) = initial_values/ C_UB-1;

gradc(:,1) = senstivity_matrix'/C_UB;




c(2) = -initial_values/C_LB - 1;


gradc(:,2) = -senstivity_matrix'/C_LB;
%


%% disp

C_LB = 0.03;
C_UB = 0.03;

[initial_values,senstivity_matrix] = read_general_sensitivity(filename,DESVAR_labels,14,'DISP');

c(3) = initial_values/C_UB -1;

gradc(:,3) = senstivity_matrix'/C_UB;




c(4) = -initial_values/C_LB-1;




gradc(:,4) = -senstivity_matrix'/C_LB;


%
c=c';
%% stress -> KS



%% 

c_stress = [];
gradc_stress =[];

C_LB = 25000;
C_UB = 25000;

[initial_values,senstivity_matrix] = read_general_sensitivity(filename,DESVAR_labels,1,'STRESS');

c_temp = initial_values/C_UB  - ones(length(initial_values),1);
c_stress=[c_stress c_temp'];

gradc_temp = senstivity_matrix'/C_UB;
gradc_stress=[gradc_stress gradc_temp];



c_temp = -initial_values/C_LB  - ones(length(initial_values),1);
c_stress=[c_stress c_temp'];

gradc_temp = -senstivity_matrix'/C_LB;
gradc_stress=[gradc_stress gradc_temp];



%%


C_LB = 25000;
C_UB = 25000;

[initial_values,senstivity_matrix] = read_general_sensitivity(filename,DESVAR_labels,2,'STRESS');

c_temp = initial_values/C_UB  - ones(length(initial_values),1);
c_stress=[c_stress c_temp'];

gradc_temp = senstivity_matrix'/C_UB;
gradc_stress=[gradc_stress gradc_temp];




c_temp = -initial_values/C_LB  - ones(length(initial_values),1);
c_stress=[c_stress c_temp'];

gradc_temp = -senstivity_matrix'/C_LB; 
gradc_stress=[gradc_stress gradc_temp];



%%


C_LB = 25000;
C_UB = 25000;

[initial_values,senstivity_matrix] = read_general_sensitivity(filename,DESVAR_labels,3,'STRESS');

c_temp = initial_values/C_UB  - ones(length(initial_values),1);
c_stress=[c_stress c_temp'];

gradc_temp = senstivity_matrix'/C_UB;
gradc_stress=[gradc_stress gradc_temp];




c_temp = -initial_values/C_LB  - ones(length(initial_values),1);
c_stress=[c_stress c_temp'];

gradc_temp = -senstivity_matrix'/C_LB;
gradc_stress=[gradc_stress gradc_temp];





%%


C_LB = 25000;
C_UB = 25000;

[initial_values,senstivity_matrix] = read_general_sensitivity(filename,DESVAR_labels,6,'STRESS');

c_temp = initial_values/C_UB  - ones(length(initial_values),1);
c_stress=[c_stress c_temp'];

gradc_temp = senstivity_matrix'/C_UB;
gradc_stress=[gradc_stress gradc_temp];




c_temp = -initial_values/C_LB  - ones(length(initial_values),1);
c_stress=[c_stress c_temp'];

gradc_temp = -senstivity_matrix'/C_LB; 
gradc_stress=[gradc_stress gradc_temp];


%% Form KS function for stress
rho = 50;
gmax = 0;

[KS,UB]=KS_opt_xAEWing3(c_stress,gmax);

c(5) = KS-UB;

g = c_stress;

nominator   = exp(rho*g)*gradc_stress';
denominator = rho*sum(exp(rho*g));

gradc(:,5)  = nominator'/denominator;
%% LAMA


C_LB = 1.5;

[initial_values,senstivity_matrix] = read_general_sensitivity(filename,DESVAR_labels,16,'LAMA');



c_temp = -initial_values/C_LB  - ones(length(initial_values),1);c=[c; c_temp];

gradc_temp = -senstivity_matrix'/C_LB; gradc=[gradc gradc_temp];




