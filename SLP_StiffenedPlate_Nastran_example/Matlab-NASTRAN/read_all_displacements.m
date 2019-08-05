function [displacements,disp_ratios] = read_all_displacements(pchfname,minodeid,maxnodeid,MAX_disp_rot)

% the displacement is punched to *.pch file



fid111=fopen(pchfname);
status=fseek(fid111,0,'eof');
EOF=ftell(fid111);
currentFPI=fseek(fid111,0,'bof');
% modeshape=zeros(maxmode,6);
stop_flag=1;

displacements=[0 0 0 0 0 0 0];  % node ID + 6 degrees of freedom
disp_ratios = [0 0 0 0 0 0];

while currentFPI<EOF && stop_flag
    
    linef06=fgetl(fid111);currentFPI=ftell(fid111);
    str=findstr(linef06,'DISPLACEMENTS');
    
    if length(linef06)>70 && isempty(str)==0
        
    linef06=fgetl(fid111);currentFPI=ftell(fid111);        
    linef06=fgetl(fid111);currentFPI=ftell(fid111);        
        
        jj=1; % The first line, jj=1, is the frequency line, set jj=2 if include frequency in the first line
        nodeid=minodeid;
        flag=1;

            
%             display(['----------------- Reading the mode shape of Mode #' num2str(mode) '----------------']);
            
            % uncomment this code if include frequency
%             eval(['modeshape' num2str(mode) '(1,1:7)=[NaturalFrequency(' num2str(mode) '),0,0,0,0,0,0];']);
            
            while currentFPI<EOF && flag && nodeid<=maxnodeid
                linef06=fgetl(fid111);currentFPI=ftell(fid111);                
                nodeid=str2num(linef06(1:16));
                if isempty(nodeid)==1 
                    flag=0;
                else
                    %                 eval(['modeshape' num2str(jj) '=mode' num2str(jj) '(5)']);
                    
                    
                   displacements(jj,1)=str2num(linef06(1:10));
                    
%                     str2num(linef06(20:72))
                    disp_xyz = str2num(linef06(20:72));
                    displacements(jj,2:4)= disp_xyz;
%                     eval(['modeshape' num2str(mode) '(' num2str(jj) ',1:4)=[nodeid,str2num(linef06(20:72))];']);
                    
                    linef06=fgetl(fid111);currentFPI=ftell(fid111);
                    
%                     str2num(linef06(20:72))
                   beta_xyz = str2num(linef06(20:72));
                   displacements(jj,5:7) = beta_xyz;
                    
                   disp_ratios(jj,1:2) = [disp_xyz(3) -disp_xyz(3)]/MAX_disp_rot(1);
                   disp_ratios(jj,3:4) = [beta_xyz(1) -beta_xyz(1)]/MAX_disp_rot(2);
                   disp_ratios(jj,5:6) = [beta_xyz(2) -beta_xyz(2)]/MAX_disp_rot(3);
%                     eval(['modeshape' num2str(mode) '(' num2str(jj) ',5:7)=str2num(linef06(20:72));']);
                    
                    jj=jj+1;
                end
            end
            
    end
end
fclose(fid111);


