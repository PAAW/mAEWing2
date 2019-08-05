function [HM,maximum_disp]=read_HM(sol144_f06)

try
    % read hinge moment for all servos from sol 144 trim analysis
    
    % v0.1
    % add three lines to read the maximum displacements because the 3 subcases
    
    HM=[];
    maximum_disp = [];
    
    fid=fopen(sol144_f06);
    
    keyword1='CONTROL';
    keyword2='SURFACE';
    keyword3='POSITION';
    % keyword4= 'AND HINGE MOMENT RESULTS'
    
    
    keyword4='MAXIMUM';
    keyword5='DISPLACEMENTS';
    
    status=fseek(fid,0,'eof');
    EOF=ftell(fid);
    currentFPI=fseek(fid,0,'bof');
    
    
    linedata=fgetl(fid);currentFPI=ftell(fid);
    
    hm_num=1;
    subcase=1;
    
    while currentFPI<EOF
        
        % read hinge moment
        if isempty(findstr(linedata,keyword1))==0 && isempty(findstr(linedata,keyword2))==0 && isempty(findstr(linedata,keyword3))==0
            
            linedata=fgetl(fid);currentFPI=ftell(fid);
            linedata=fgetl(fid);currentFPI=ftell(fid);
            linedata=fgetl(fid);currentFPI=ftell(fid);
            linedata=fgetl(fid);currentFPI=ftell(fid);
            linedata=fgetl(fid);currentFPI=ftell(fid);
            
            
            linedata=fgetl(fid);currentFPI=ftell(fid);
            
            while length(linedata)==127 %&& abs(str2num(linedata(30:43))+10/180*pi)<1e-4
                
                
                HM(hm_num,1)=str2num(linedata(31:43));
                HM(hm_num,2)=str2num(linedata(46:59));
                HM(hm_num,3)=str2num(linedata(64:75));
                HM(hm_num,4)=str2num(linedata(82:94));
                HM(hm_num,5)=str2num(linedata(97:110));
                HM(hm_num,6)=str2num(linedata(114:126));
                
                
                
                linedata=fgetl(fid);currentFPI=ftell(fid);
                
                hm_num=hm_num+1;
            end
            
            
            
            % read maximum displacement
            
            
        elseif isempty(findstr(linedata,keyword4))==0 && isempty(findstr(linedata,keyword5))==0
            
            linedata=fgetl(fid);currentFPI=ftell(fid);
            linedata=fgetl(fid);currentFPI=ftell(fid);
            
            linedata=fgetl(fid);currentFPI=ftell(fid);
            linedata_num=str2num(linedata);
            maximum_disp(1)=linedata_num(5);
            
            %         linedata=fgetl(fid);currentFPI=ftell(fid);
            %         linedata_num=str2num(linedata);
            %         maximum_disp(2)=linedata_num(5);
            %
            %         linedata=fgetl(fid);currentFPI=ftell(fid);
            %         linedata_num=str2num(linedata);
            %         maximum_disp(3)=linedata_num(5);
            %         subcase=subcase+1;
        end
        
        % read again
        linedata=fgetl(fid);currentFPI=ftell(fid);
    end
    
    
    if isempty(HM)==1
        
        HM=[10000 10000 10000 10000 10000 10000];
        
    end
    
    fclose(fid);
    
catch
    
    HM=[10000 10000 10000 10000 10000 10000];
    maximum_disp=1000;
    
end
