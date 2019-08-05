function [f,gradf] = fitness(x0)


% save x0-fitness.txt x0 -ascii -append

%% WRITE DESVAR.BDF


updated_DESVAR(x0);


%% RUN NASTRAN SOL 200

design_folder = 'C:\Users\weizhao\Documents\PAAW\N+3\MATLAB_NASTRAN_Opt\verification_stiffened_plate';
main_folder = 'C:\Users\weizhao\Documents\PAAW\N+3\MATLAB_NASTRAN_Opt\sensitivity_extract';

nastran = 'C:\MSC.Software\MSC_Nastran\20170\bin\nast20170.exe ';


sol200_bdf_fname='dsoug4_only_sens';
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
DESVAR_labels = {'T-PLATE'
    'HATDIM2'};
filename  = [design_folder filesep sol200_bdf_fname '.f06'];

[initial_values,senstivity_matrix] = read_general_sensitivity(filename,DESVAR_labels,15,'WEIGHT');

f = initial_values;

gradf = senstivity_matrix;


save weight_iteration.txt f -ascii -append