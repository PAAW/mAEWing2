function [FI,ST]=readFI_opt_v2(sol144_f06)

% read failure index


%'  20600095   STRAIN      1            0.0623   -12';
%'  10200072   TSAI-WU     1            0.0014      ';

try
    FI =[];
    ST=[];
    fid111=fopen(sol144_f06);
    status=fseek(fid111,0,'eof');
    EOF=ftell(fid111);
    currentFPI=fseek(fid111,0,'bof');
    
    FT1 = 'TSAI-WU';
    FT2 = 'STRAIN';
    
    TW_FI_num = 1;
    ST_FI_num = 1;
    
    
    while currentFPI<EOF
        
        linef06=fgetl(fid111);currentFPI=ftell(fid111);
        
        if findstr(linef06,FT1)
            
            TW_flag = 1;
            
            while TW_flag
                
                linef06=fgetl(fid111);currentFPI=ftell(fid111);
                
                if length(linef06)==124
                    
                    if isempty(linef06(100:end))==0
                        %    linef06
                        max_elem_FI=str2num(linef06(101:110));
                        
                        if isempty( max_elem_FI)==0 
                            
                            FI(TW_FI_num)=max_elem_FI(1);
                            
                            TW_FI_num=TW_FI_num+1;
                             TW_flag = 0;
                        end
                    end
                   
                end

            end
            
            
        elseif findstr(linef06,FT2)
            
            Max_STRAIN_flag = 1;
            
            while Max_STRAIN_flag
                
               
                linef06=fgetl(fid111);
                currentFPI=ftell(fid111);
                
                if length(linef06)==124
                    
                    if isempty(linef06(100:end))==0
                        %    linef06
                        max_elem_ST=str2num(linef06(101:110));
%                          max_elem_ST
                        if isempty( max_elem_ST)==0
                            
                            ST(ST_FI_num)=max_elem_ST(1);
                            
                            ST_FI_num=ST_FI_num+1;
                             Max_STRAIN_flag = 0;
                        end
                    end
                   
                end
                
                
            end
            
            
        end
    end
    
    fclose(fid111);
    

catch
    
    ST = 1E4;
    FI=1E4;
end
