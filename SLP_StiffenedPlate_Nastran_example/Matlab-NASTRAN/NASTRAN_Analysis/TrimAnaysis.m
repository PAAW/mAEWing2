function [flap_deflections,flap_hingemoments,AoA,MASS,max_disp,FI]=...
    TrimAnaysis(design_folder,design_file_name,input_value,flagER,PATH)

% This is for trim analysis of mAEWing2 with discrete flaps


% flag=1: Elastic trim analysis
% flag=0: Rigid trim analysis



% input_value: Trailing-edge flap deflections + Leading-edge flap
% deflections Delta_TF + Delta_LF


TE_flaps=input_value(1:length(input_value)-1);

owleadL=input_value(length(input_value));

% copy
sol144_bdf_base='Sol144_mAEWing2_1G.bdf';

copyfile([PATH.base_bdf_folder '\' sol144_bdf_base],[design_folder '\sol144_' design_file_name '.bdf' ]);


fid102=fopen([design_folder '\sol144_' design_file_name '.bdf' ],'a+');
% open with append


TE_flaps_num=length(input_value)-1;


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
    line_bdf=['TRIM ,   101   ,  0.1    , .1046 ,  URDD3,	1.0    , owleadL, ' num2str8(owleadL) ', 0.0,  +TR1'];
elseif flagER==1
    line_bdf=['TRIM ,   101   ,  0.1    , .1046 ,  URDD3,	1.0    , owleadL ,' num2str8(owleadL) ',  ,  +TR1'];
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
call_NASTRAN_sol144_discrete=[PATH.nastran ' ' design_folder '\sol144_' design_file_name '.bdf' ];
system(call_NASTRAN_sol144_discrete);

cd(PATH.mainprogram_folder);


% read flap deflections, AoA, leading flap deflection, and hinge moments

sol144_f06=[design_folder '\sol144_' design_file_name '.f06' ];

%% READ Trim analysis results in hinge moments and flap deflections

[HM,max_disp]=read_HM(sol144_f06);

flap_deflections=HM(:,2);

flap_hingemoments=HM(:,5);

% READ FI
FI=readFI(sol144_f06,1e-3);
% read trimmed AoA

fid_fort54=fopen([design_folder '\fort.54']);
linedata=fgetl(fid_fort54);

% case 1
linedata=fgetl(fid_fort54);
trimmed_data=str2num(linedata);
AoA(1)=trimmed_data(9);
% % case 2
% linedata=fgetl(fid_fort54);
% trimmed_data=str2num(linedata);
% AoA(2)=trimmed_data(9)/pi*180;
% % case 3
% linedata=fgetl(fid_fort54);
% trimmed_data=str2num(linedata);
% AoA(3)=trimmed_data(9)/pi*180;


MASS.WEIGHT=trimmed_data(21);
MASS.XCG=trimmed_data(22);


fclose(fid_fort54);




%% READ Failure index


