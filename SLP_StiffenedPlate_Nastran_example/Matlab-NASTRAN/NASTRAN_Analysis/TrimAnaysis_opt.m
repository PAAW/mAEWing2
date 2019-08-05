function [flap_deflections,flap_hingemoments,AoA,WEIGHT,max_disp,FI]=...
    TrimAnaysis_opt(design_folder,...
    TE_flaps,LE_flap,flagER,base_bdf_folder,nastran, main_program_folder,...
    URDD3,dynamicpressure,max_Hingemoment,AoA_max,BF_max,CL_max,...
    ref_area,ref_chord,TE_flap_max,LE_flap_max)

% This is for trim analysis of mAEWing2 with discrete flaps


% flag=1: Elastic trim analysis
% flag=0: Rigid trim analysis


% input_value: Trailing-edge flap deflections + Leading-edge flap
% deflections Delta_TF + Delta_LF


% TE_flaps=TE_LE_flap(1:length(TE_LE_flap)-1);

% LE_flap=TE_LE_flap(length(TE_LE_flap));

% copy


%% DETERMINE THE OPTIMAL FLAP

num_of_flaps = length(TE_flaps)+length(LE_flap);

x_0 = [TE_flaps,LE_flap, URDD3]';

Optimal_Flap = FlapOptimization(x_0,dynamicpressure,...
    design_folder,main_program_folder,base_bdf_folder,nastran,...
    max_Hingemoment,ref_area,ref_chord,...
    AoA_max,BF_max,CL_max,URDD3,...
    TE_flap_max,LE_flap_max);

%%

disp('====== R U N            O P T I M A L        F L A P =======');
TE_flaps = Optimal_Flap(1:num_of_flaps-1);
LE_flap = Optimal_Flap(num_of_flaps);
sol144_bdf_base='sol144.bdf';

copyfile([base_bdf_folder filesep sol144_bdf_base],[design_folder filesep 'sol144.bdf' ]);


fid102=fopen([design_folder filesep 'sol144.bdf' ],'a+');
% open with append


TE_flaps_num=length(TE_flaps);


for TE_flap_num=1:TE_flaps_num
    
    
    line_bdf=['AELINK, ALWAYS, owTER' num2str(TE_flap_num) ',  owTEL'  num2str(TE_flaps_num+1-TE_flap_num) ', -1.0'];
    
    fprintf(fid102,'%s\n',line_bdf);
end


line_bdf='$====================================================';
fprintf(fid102,'%s\n',line_bdf);

line_bdf='$                   TRIM ANALYSIS';
fprintf(fid102,'%s\n',line_bdf);

line_bdf='$====================================================';
fprintf(fid102,'%s\n',line_bdf);

line_bdf='$1......2.......3.......4.......5.......6.......7.......8.......9.......0';
fprintf(fid102,'%s\n',line_bdf);



% check how many lines in TRIM card based on number of trailing-edge flaps
%% TRIM=101; 1-G LOAD
fprintf(fid102,'%s\n','$ ---- 1G LOAD ----');
if flagER==0
    line_bdf=['TRIM ,   101   ,  0.1    , ' num2str8(dynamicpressure) ' ,  URDD3,' num2str8(URDD3) ', owleadL, ' num2str8(LE_flap) ', 0.0,  +TR1'];
elseif flagER==1
    line_bdf=['TRIM ,   101   ,  0.1    ,' num2str8(dynamicpressure) ' ,  URDD3,' num2str8(URDD3) ', owleadL ,' num2str8(LE_flap) ',  ,  +TR1'];
end
fprintf(fid102,'%s\n',line_bdf);

line_bdf=['+TR1 ,   PITCH  , 0.00,URDD5 ,0.0 , owTEL1 ,  ' num2str8(TE_flaps(1)) '   ,  owTEL2 , ' num2str8(TE_flaps(2)) ' ,   +TR2'];
fprintf(fid102,'%s\n',line_bdf);

line_bdf=['+TR2 , owTEL3 ,' num2str8(TE_flaps(3)) ', owTEL4 ,' num2str8(TE_flaps(4))];

fprintf(fid102,'%s\n',line_bdf);

% %% TRIM 102; -1G LOAD
% fprintf(fid102,'%s\n','$ ---- -1G LOAD ----');
% if flagER==0
%     line_bdf=['TRIM ,   102   ,  0.1    , .1046 ,  URDD3,	-1.0    , owleadL, ' num2str8(owleadL) ', 0.0,  +TR1'];
% elseif flagER==1
%     line_bdf=['TRIM ,   102   ,  0.1    , .1046 ,  URDD3,	-1.0    , owleadL ,' num2str8(owleadL) ',  ,  +TR1'];
% end
% fprintf(fid102,'%s\n',line_bdf);
% 
% line_bdf=['+TR1 ,   PITCH  , 0.00,URDD5 ,0.0 , owTEL1 ,  ' num2str8(TE_flaps(1)) '   ,  owTEL2 , ' num2str8(TE_flaps(2)) ' ,   +TR2'];
% fprintf(fid102,'%s\n',line_bdf);
% 
% line_bdf=['+TR2 , owTEL3 ,' num2str8(TE_flaps(3)) ', owTEL4 ,' num2str8(TE_flaps(4))];
% 
% fprintf(fid102,'%s\n',line_bdf);
% 
% 
% 
% %% TRIM 103; +3.8G load
% 
% fprintf(fid102,'%s\n','$ ---- +3.8G LOAD ----');
% if flagER==0
%     line_bdf=['TRIM ,   103   ,  0.1    , .1046 ,  URDD3,	3.8    , owleadL, ' num2str8(owleadL) ', 0.0,  +TR1'];
% elseif flagER==1
%     line_bdf=['TRIM ,   103   ,  0.1    , .1046 ,  URDD3,	3.8    , owleadL ,' num2str8(owleadL) ',  ,  +TR1'];
% end
% fprintf(fid102,'%s\n',line_bdf);
% 
% line_bdf=['+TR1 ,   PITCH  , 0.00,URDD5 ,0.0 , owTEL1 ,  ' num2str8(TE_flaps(1)) '   ,  owTEL2 , ' num2str8(TE_flaps(2)) ' ,   +TR2'];
% fprintf(fid102,'%s\n',line_bdf);
% 
% line_bdf=['+TR2 , owTEL3 ,' num2str8(TE_flaps(3)) ', owTEL4 ,' num2str8(TE_flaps(4))];
% 
% fprintf(fid102,'%s\n',line_bdf);


%%


line_bdf='$ output lift distribution';
fprintf(fid102,'%s\n',line_bdf);


line_bdf='MONCNCM ALL     ALL STRIPS';fprintf(fid102,'%s\n',line_bdf);
line_bdf='                all';fprintf(fid102,'%s\n',line_bdf);
line_bdf='ENDDATA d2dbc08c';fprintf(fid102,'%s\n',line_bdf);

fclose(fid102);


% RUN NASTRAN

%% RUN SOL144 AT 1-G LEVEL FLIGHT
% cd([PATH.design_folder '\' design_file_name]);
cd(design_folder);
call_NASTRAN_sol144_discrete=[nastran ' ' design_folder filesep 'sol144.bdf' '  memory=estimate'];
system(call_NASTRAN_sol144_discrete);
% !del *.xdb

% check program end or not
if isunix
    fname = [design_folder filesep 'sol144'];
    program_nastran_determine(fname)
end
%
cd(main_program_folder);


% read flap deflections, AoA, leading flap deflection, and hinge moments

sol144_f06=[design_folder filesep 'sol144.f06' ];

%% READ Trim analysis results in hinge moments and flap deflections

[HM,max_disp]=read_HM(sol144_f06);

flap_deflections=HM(:,2);

flap_hingemoments=HM(:,5);

% READ FI
FI=readFI_opt(sol144_f06,1e-3);
% read trimmed AoA

try
    fid_fort54=fopen([design_folder filesep 'fort.54']);
    linedata=fgetl(fid_fort54);
    
    % case 1
    linedata=fgetl(fid_fort54);
    trimmed_data=str2num(linedata);
    AoA(1)=trimmed_data(9);
    WEIGHT=trimmed_data(21);
    
catch
    AoA=pi;
    WEIGHT=1E4;
end
% % case 2
% linedata=fgetl(fid_fort54);
% trimmed_data=str2num(linedata);
% AoA(2)=trimmed_data(9)/pi*180;
% % case 3
% linedata=fgetl(fid_fort54);
% trimmed_data=str2num(linedata);
% AoA(3)=trimmed_data(9)/pi*180;


% MASS.XCG=trimmed_data(22);

fclose(fid_fort54);




%% READ Failure index


