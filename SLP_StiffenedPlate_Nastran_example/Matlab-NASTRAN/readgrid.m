function GridCord=readgrid(filename,CORD)
% ************************************************
% Read grid ID and coordinate from *f06 file with ECHO=SORT
% Grid = [nodelabel, xcord, ycord, zcord].
% read local coordinate system and output coordinate defined in global
% coordinate system, Sep,19,2016

fid111=fopen(filename);
status=fseek(fid111,0,'eof');
EOF=ftell(fid111);
currentFPI=fseek(fid111,0,'bof');
jj=1;
while currentFPI<EOF
    linef06=fgetl(fid111);currentFPI=ftell(fid111);
    str=findstr(linef06,'GRID');
    if length(linef06)>70 && isempty(str)==0 
%         disp('---- find the correct line ----');
        if length(linef06)>77
            gridata=linef06(39:77);
        else
            gridata=linef06(39:end);
        end
        
        nodenum=gridata(1:8);
        cordID = gridata(9:16);
        
        
        if strcmp(cordID,'        ') || cordID(1) =='0'
            
            deltaX=0;
            deltaY=0;
            deltaZ=0;
            
                    
        gridata1=gridata(17:24);
        gridata2=gridata(25:32);
        gridata3=gridata(33:end);
        gridcord(jj,:)=[str2num(nodenum),str2num(gridata1)+deltaX,str2num(gridata2)+deltaY,str2num(gridata3)+deltaZ];
            
        else 
            
            cord_num=find(CORD(:,1)==str2num(cordID));
            
            deltaX=CORD(cord_num,2);
            deltaY=CORD(cord_num,3);
            deltaZ=CORD(cord_num,4);
            
            Delta=[deltaX deltaY deltaZ];
            
            Theta = CORD(cord_num,5);
                    
        gridata1=gridata(17:24);
        gridata2=gridata(25:32);
        gridata3=gridata(33:end);
        str2num(nodenum);
        OriginalCord=[str2num(gridata1),str2num(gridata2),str2num(gridata3)];
        
        Coordinates = Transform_RZ(Delta,Theta,OriginalCord);
        
        gridcord(jj,:)=[str2num(nodenum),Coordinates'];
        
        
        end

%         gridcordLabel(jj,:)=str2num(nodenum)
        jj=jj+1;
    end
end
GridCord=gridcord;
fclose(fid111);