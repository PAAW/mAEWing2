function program_nastran_determine(fname)


flagENDlog=0;

while flagENDlog==0
    
    pause(10); % allow nastran to generate *.log file
    
    tic
    
    
    
    if exist([fname '.log'],'file')==2
        
        flag_finished = 0;
        flag_time_limit = 1;
        
        
        while flag_finished==0 && flag_time_limit
            
            fid100=fopen([fname '.log']);
            status=fseek(fid100,0,'eof');
            EOF=ftell(fid100);
            currentFPI=fseek(fid100,0,'bof');
            while currentFPI<EOF
                str=fgetl(fid100);
                currentFPI=ftell(fid100);
                if length(str)>45
                    if strcmp(str(1:20),'MSC Nastran finished')==1
                        flag_finished=1;
                        
                        disp('--- Successfully finish program running ----')
                        %                         break
                    end
                    
                end
            end
            
            fclose(fid100);
            
            %% Time
            elapsed_time=toc;
            
            if elapsed_time>3000 % sol145 takes around 1000/60
                
                % no log file generated, stop this program
                disp('--- Exceed time limit ----')
                flag_time_limit=0;
            end
            
            
            % %             flag_time_limit
            % %             flag_finished
            
        end
        
    end
    
    flagENDlog=1;
    
end
