% read lift distribution
function [YS,ZS,XS,REFC,REFS,RigidLift,ElasticLift,ElasticMoment]=read_trim_results_opt_v2(filename)
% change 9 - to 8 in YS read to captuer the negative vlaue of span location
% read TRIM SOL144 Results for TrimOpt

% try


string1 = 'STRP';
string2 = 'SUBCASE';
fid111=fopen(filename);
status=fseek(fid111,0,'eof');
EOF=ftell(fid111);
currentFPI=fseek(fid111,0,'bof');
% subcase=1;
% flag=1;

while currentFPI<EOF
    
    linef06=fgetl(fid111);currentFPI=ftell(fid111);
    
    
    subcase_id = findstr(linef06,string2);
    
    if isempty(subcase_id)==0

        subcase = str2num(linef06(subcase_id+7:end));
    end
    
    strp=findstr(linef06,string1);
    
    strip_num=0;
    
    if isempty(strp)==0
        
        %             subcase=subcase+1;
        flag=1;
        
        linef06=fgetl(fid111);currentFPI=ftell(fid111);
        linef06=fgetl(fid111);currentFPI=ftell(fid111);
        
        strip_number_old=strip_num;
        
        %         str2num(linef06(1:5))
        
        while flag
            
            strip_number_new=str2num(linef06(1:5));
            
            YS(strip_number_new,subcase)=str2num(linef06(8:18));
            
            ZS(strip_number_new,subcase)=str2num(linef06(20:30));
            
            XS(strip_number_new,subcase)=str2num(linef06(31:43));
            
            REFC(strip_number_new,subcase)=str2num(linef06(44:55));
            
            REFS(strip_number_new,subcase)=str2num(linef06(57:68));
            
            RigidLift(strip_number_new,subcase)=str2num(linef06(70:81));
            
            ElasticLift(strip_number_new,subcase)=str2num(linef06(83:94));
            
            ElasticMoment(strip_number_new,subcase)=str2num(linef06(109:end));
            
            linef06=fgetl(fid111);currentFPI=ftell(fid111);
            
            strip_number_old=strip_number_new;
            
            if length(linef06)>=100 && str2num(linef06(1:5))==strip_number_old+1
                flag=1;
            else
                
                flag=0;
            end
        end
    end
    
    
    %     linef06=fgetl(fid111);currentFPI=ftell(fid111);
    
end

fclose(fid111);
% catch
%
%     YS=[];
%     REFC=[];
%     REFS=[];
%     RigidLift=[];
%     ElasticLift=[];


