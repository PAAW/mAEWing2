function [max_strain,max_disp]=xAEWing3_sol144_trim(base_bdf_folder,design_folder,...
    flaps,load_factors,dynamic_pressure,...
    sol144_bdf_base_fname,sol144_bdf_fname,...
    ref_chord,ref_area,...
    TE_flaps_num,Mach,...
    nastran,main_program_folder)

% run Trim analysis for optimal flap rotations in 1G, 3.8G and lift
% maximization case



subcases_num = size(flaps,1);

max_strain = zeros(subcases_num,1);

max_disp = zeros(subcases_num,1);


%%
copyfile([base_bdf_folder filesep sol144_bdf_base_fname '.bdf'],[design_folder filesep sol144_bdf_fname '.bdf' ]);
fid102=fopen([design_folder filesep sol144_bdf_fname '.bdf'],'a+');



for case_num = 1:subcases_num
    
    case_line=['SUBCASE ' num2str(case_num)];
    fprintf(fid102,'%s\n',case_line);
    case_line=['   TRIM =  ' num2str(case_num+100)];
    fprintf(fid102,'%s\n',case_line);
    
     fprintf(fid102,'%s\n','   TRIMF = ALL');
    
end


case_line='BEGIN BULK';
fprintf(fid102,'%s\n',case_line);
case_line='INCLUDE  ''xAEWing3_structure_aero.dat''';
fprintf(fid102,'%s\n',case_line);

% AERO

ref_span = ref_area/ref_chord;
aeros_line = ['AEROS   0       0       ' num2str8(ref_chord) num2str8(ref_span) num2str8(ref_area)];
fprintf(fid102,'%s\n',aeros_line);

% LINE
line_bdf='$====================================================';
fprintf(fid102,'%s\n',line_bdf);
line_bdf='$                   TRIM ANALYSIS';
fprintf(fid102,'%s\n',line_bdf);
line_bdf='$====================================================';
fprintf(fid102,'%s\n',line_bdf);
line_bdf='$1......2.......3.......4.......5.......6.......7.......8.......9.......0';
fprintf(fid102,'%s\n',line_bdf);

%% AELINK

for TE_flap_num=1:TE_flaps_num
    
    line_bdf=['AELINK, ALWAYS, owTER' num2str(TE_flap_num) ',  owTEL'  num2str(TE_flaps_num+1-TE_flap_num) ', -1.0'];
    
    fprintf(fid102,'%s\n',line_bdf);
end




%% create sol 144 with multiple load cases
flagER = 1;
for x_num = 1:subcases_num
    
    
    %% ================
    URDD3 = load_factors(x_num);
    xflap = flaps(x_num,:);
    
    %% TRIM CARD
    TrimID = 100 + x_num ;
    write_TRIM_many_flaps_mAEWing2_opt_v2(dynamic_pressure(x_num),URDD3,xflap,TE_flaps_num,fid102,flagER,Mach,TrimID)
    
end

%% output lift distribution
% fprintf(fid102,'%s\n','$$ ==== OUTPUT LIFT DISTRIBUTION ====');
% line_bdf='MONCNCM ALL     ALL STRIPS';fprintf(fid102,'%s\n',line_bdf);
% line_bdf='                all';fprintf(fid102,'%s\n',line_bdf);
line_bdf='ENDDATA d2dbc08c';fprintf(fid102,'%s\n',line_bdf);
% 
fclose(fid102);

%% RUN nastran

cd(design_folder);
call_NASTRAN_sol144_discrete=[nastran ' ' design_folder filesep sol144_bdf_fname '.bdf' '  memory=2g'];
system(call_NASTRAN_sol144_discrete);
!del *.xdb

% check program end or not
if isunix
    fname = [design_folder filesep sol144_bdf_fname];
    program_nastran_determine(fname)
end
%
cd(main_program_folder);

%% Read maximum displacement and maximum failure index
sol144_f06 = [design_folder filesep sol144_bdf_fname '.f06'];





