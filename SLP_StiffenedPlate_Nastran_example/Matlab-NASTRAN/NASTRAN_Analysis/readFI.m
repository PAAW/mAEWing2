function FI=readFI(sol144_f06,Tol)

% read failure index


try
    
    fid111=fopen(sol144_f06);
    
    
    status=fseek(fid111,0,'eof');
    EOF=ftell(fid111);
    currentFPI=fseek(fid111,0,'bof');
    
    
    STR='TSAI-WU';
    
    FI_max_num=1;
    
    while currentFPI<EOF
        
        linef06=fgetl(fid111);currentFPI=ftell(fid111);
        
        
        
        if length(linef06)==124
            
            if isempty(str2num(linef06(100:end)))==0
                
                
                %             linef06
                max_elem_FI=str2num(linef06);
                
                if max_elem_FI(1)>Tol
                    
                    FI(FI_max_num)=max_elem_FI(1);
                    
                    
                    
                    
                    FI_max_num=FI_max_num+1;
                    
                end
                
            end
            
            
            
            %     FI_flap=strfind(linef06,STR);
            %
            %     if isempty(FI_flap)==0
            %
            %
            %         linef06=fgetl(fid111);currentFPI=ftell(fid111);
            %         linef06=fgetl(fid111);currentFPI=ftell(fid111);
            %         linef06=fgetl(fid111);currentFPI=ftell(fid111);
            %         linef06=fgetl(fid111);currentFPI=ftell(fid111);
            %         linef06=fgetl(fid111)
            %         currentFPI=ftell(fid111);
            %
            %
            %
            %         max_elem_FI=str2num(linef06);
            %         FI(FI_max_num)=max_elem_FI(1);
            %
            %         FI_max_num=FI_max_num+1;
            %     end
            
            
            
            
        end
        
    end
    
    fclose(fid111);
    
catch
    
    
    FI=1E4;
end
