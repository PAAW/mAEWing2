% read stress senstivity% stress components
function [stress0,stress_sens] = read_stress_components_sensitivity(filename,DESVAR_labels,DRESP1_ID)

%% Only use DRESP1 here, the senstivity of DRESP2 can be obtained from DRESP1

% DRESP1_ID = 2; % if specify properties, there would be many elements.

% filename = 'C:\Users\weizhao\Documents\PAAW\N+3\MATLAB_NASTRAN_Opt\verification_stiffened_plate\dsoug4_only_sens.f06';


fid100=fopen(filename);

%%
% DESVAR_labels = {'T-PLATE'
%     'HATDIM2'};


% total_subcases = 2;

stress0 = [];
stress_sens = [];

%% ===== STRESS =======

keywords2 = 'DRESP1 ID';

keywords3 = 'RESPONSE TYPE';

status=fseek(fid100,0,'eof');
EOF=ftell(fid100);
currentFPI=fseek(fid100,0,'bof');


while currentFPI<EOF
    
    linebdf=fgetl(fid100);currentFPI=ftell(fid100);
    
    if length(linebdf)>120 && isempty(findstr(linebdf,keywords2)) == 0  && isempty(findstr(linebdf,keywords3)) == 0  && isempty(findstr(linebdf,'STRESS')) == 0
        
        dresp1_ID =  str2num(linebdf(15:24));
        
        dresp1_index = find(dresp1_ID==DRESP1_ID);
        
        while length(linebdf)>1 && isempty(dresp1_index )==0
            
            linebdf=fgetl(fid100);  currentFPI=ftell(fid100);
            
            if length(linebdf)>45 && isempty(str2num(linebdf(37:45))) == 0
                
                
                
                % subcase
                if isempty(str2num(linebdf(5:11))) == 0
                    
                    subcase = str2num(linebdf(5:11));
                    
                    %                         elseif isempty(string_temp) == 0
                end
                
                % initial value
                if isempty(str2num(linebdf(14:24))) == 0
                    
                    
                    
                    
                    stress0 = [stress0 ;str2num(linebdf(14:24))];
                    
                    %                         elseif isempty(string_temp) == 0
                end
                
                
                
                %
                
                string_temp = linebdf(38:end);
                
                for desvar_label  = 1:size(DESVAR_labels,1)
                    
                    desvar_temp = DESVAR_labels{desvar_label};
                    
                    desvar_exit = findstr( desvar_temp,string_temp);
                    
                    if isempty(desvar_exit) == 0
                        
                        stress_sens_temp (desvar_label) = str2num(string_temp(desvar_exit +10: desvar_exit+20));
                        
                    end
                    
                end
                
                
                stress_sens =[ stress_sens;    stress_sens_temp];
                
                
            end
            
        end
        
        
    end
    
end


fclose('all');


