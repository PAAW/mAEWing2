%% READ NASTRAN SOL 200 Sensitivity Matrix

%% Only use DRESP1 here, the senstivity of DRESP2 can be obtained from DRESP1
filename = 'C:\Users\weizhao\Documents\PAAW\N+3\MATLAB_NASTRAN_Opt\verification_stiffened_plate\dsoug4_only_sens.f06';


fid100=fopen(filename);


stress_sens = [];

disp_sens = [];

weight_sens =[];

%%
DESVAR_labels = {'T-PLATE'
    'HATDIM2'};


weight_sens = zeros(1,size(DESVAR_labels,1));

total_subcases = 2;





%% === Disp ==== need to add component

DISP_GRID_ID = [10302
    10203];

disp_constraints_number = size(DISP_GRID_ID,1);

disp_sen = zeros(disp_constraints_number,size(DESVAR_labels,1),total_subcases);
disp_values = zeros(disp_constraints_number,total_subcases);


[initial_values,senstivity_matrix] = read_general_sensitivity(filename,DESVAR_labels,13,'DISP');


[initial_values,senstivity_matrix] = read_general_sensitivity(filename,DESVAR_labels,14,'DISP');

%% ===== STRESS =======

DRESP1_ID = 1

[stress0,stress_sens] = read_stress_components_sensitivity(filename,DESVAR_labels,DRESP1_ID);


DRESP1_ID = 2

[stress0,stress_sens] = read_stress_components_sensitivity(filename,DESVAR_labels,DRESP1_ID);


DRESP1_ID = 3

[stress0,stress_sens] = read_stress_components_sensitivity(filename,DESVAR_labels,DRESP1_ID);


DRESP1_ID = 6

[stress0,stress_sens] = read_stress_components_sensitivity(filename,DESVAR_labels,DRESP1_ID);

%%
if fid100>0
    
    %     keywords1 = '*    R E S P O N S E    S E N S I T I V I T Y    C O E F F I C I E N T S   *';
    keywords2 = 'DRESP1 ID';
    
    keywords3 = 'RESPONSE TYPE';
    
    status=fseek(fid100,0,'eof');
    EOF=ftell(fid100);
    currentFPI=fseek(fid100,0,'bof');
    
    
    while currentFPI<EOF
        
        linebdf=fgetl(fid100);currentFPI=ftell(fid100);
        
        if length(linebdf)>120 && isempty(findstr(linebdf,keywords2)) == 0  && isempty(findstr(linebdf,keywords3)) == 0
              %% ======== Read sensitivity for weight ==========
            if isempty(findstr(linebdf,'WEIGHT')) == 0 % No subcase
                
                % to read many design variables
                
                while length(linebdf)>1
                    
                    linebdf=fgetl(fid100);  currentFPI=ftell(fid100);
                    
                    
                    
                    if length(linebdf)>45 && isempty(str2num(linebdf(37:45))) == 0
                        
                        if isempty(str2num(linebdf(14:24))) == 0  % initial value
                            
                            weight_present_value = str2num(linebdf(14:24)) % delievered to f
                            
                            %                         elseif isempty(string_temp) == 0
                        end
                        
                        string_temp = linebdf(38:end);
                        for desvar_label  = 1:size(DESVAR_labels,1)
                            
                            desvar_temp = DESVAR_labels{desvar_label};
                            
                            desvar_exit = findstr( desvar_temp,string_temp);
                            
                            if isempty(desvar_exit) == 0
                                
                                weight_sens(desvar_label) = str2num(string_temp(desvar_exit +10: desvar_exit+20));
                                
                            end
                            
                        end
                        
                    end
                    
                    
                end
                
                
                
                %% ======== Read sensitivity for Displacement ==========
            elseif isempty(findstr(linebdf,'DISP')) == 0
                
                disp_constrain_ID =  find(str2num(linebdf(74:88)) == DISP_GRID_ID);
                
                
                
                while length(linebdf)>1
                    
                    linebdf=fgetl(fid100);  currentFPI=ftell(fid100);
                    
                    if length(linebdf)>45 && isempty(str2num(linebdf(37:45))) == 0
                        
                        if isempty(str2num(linebdf(5:11))) == 0  % initial value
                            
                            subcase = str2num(linebdf(5:11));
                            
                            %                         elseif isempty(string_temp) == 0
                        end
                        
                        
                        if isempty(str2num(linebdf(14:24))) == 0  % initial value
                            
                            disp_values(disp_constrain_ID,  subcase)= str2num(linebdf(14:24));
                            
                            %                         elseif isempty(string_temp) == 0
                        end
                        
                        
                        
                        string_temp = linebdf(38:end);
                        
                        for desvar_label  = 1:size(DESVAR_labels,1)
                            
                            desvar_temp = DESVAR_labels{desvar_label};
                            
                            desvar_exit = findstr( desvar_temp,string_temp);
                            
                            if isempty(desvar_exit) == 0
                                
                                disp_sens(disp_constrain_ID,desvar_label,subcase) = str2num(string_temp(desvar_exit +10: desvar_exit+20));
                                
                            end
                            
                        end
                        
                    end
                    
                end
                
                  %% ======== Read sensitivity for stress ==========
                
            elseif isempty(findstr(linebdf,'STRESS')) == 0
                
                                disp_constrain_ID =  find(str2num(linebdf(74:88)) == DISP_GRID_ID);
                
                
                
                while length(linebdf)>1
                    
                    linebdf=fgetl(fid100);  currentFPI=ftell(fid100);
                    
                    if length(linebdf)>45 && isempty(str2num(linebdf(37:45))) == 0
                        
                        if isempty(str2num(linebdf(5:11))) == 0  % subcase
                            
                            subcase = str2num(linebdf(5:11));
                            
                            %                         elseif isempty(string_temp) == 0
                        end
                        
                        
                        if isempty(str2num(linebdf(14:24))) == 0  % initial value
                            
                            disp_values(disp_constrain_ID,  subcase)= str2num(linebdf(14:24));
                            
                            %                         elseif isempty(string_temp) == 0
                        end
                        
                        
                        
                        string_temp = linebdf(38:end);
                        
                        for desvar_label  = 1:size(DESVAR_labels,1)
                            
                            desvar_temp = DESVAR_labels{desvar_label};
                            
                            desvar_exit = findstr( desvar_temp,string_temp);
                            
                            if isempty(desvar_exit) == 0
                                
                                disp_sens(disp_constrain_ID,desvar_label,subcase) = str2num(string_temp(desvar_exit +10: desvar_exit+20));
                                
                            end
                            
                        end
                        
                    end
                    
                end
                
                
                
                
            end
            
            
        end
        
    end
    
    
    
end




fclose('all');


disp_sens


disp_values





