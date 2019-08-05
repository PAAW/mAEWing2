% read lift distribution
function [YS_updated,ZS_updated,XS_updated,REFC_updated,REFS_updated,RigidLift_updated,ElasticLift_updated,ElasticMoment_updated]=read_trim_results(filename)
% change 9 - to 8 in YS read to captuer the negative vlaue of span location
% read TRIM SOL144 Results for TrimOpt

% try


string1='STRP';

fid111=fopen(filename);
status=fseek(fid111,0,'eof');
EOF=ftell(fid111);
currentFPI=fseek(fid111,0,'bof');
subcase=1;
% flag=1;

while currentFPI<EOF
    
    linef06=fgetl(fid111);currentFPI=ftell(fid111);
    
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
            
            YS_updated(strip_number_new,subcase)=str2num(linef06(8:18));
            
            ZS_updated(strip_number_new,subcase)=str2num(linef06(20:30));
            
            XS_updated(strip_number_new,subcase)=str2num(linef06(31:42));
            
            REFC_updated(strip_number_new,subcase)=str2num(linef06(44:55));
            
            REFS_updated(strip_number_new,subcase)=str2num(linef06(57:68));
            
            RigidLift_updated(strip_number_new,subcase)=str2num(linef06(70:81));
            
            ElasticLift_updated(strip_number_new,subcase)=str2num(linef06(83:94));
            
            ElasticMoment_updated(strip_number_new,subcase)=str2num(linef06(109:end));
            
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

%%


% [YS_updated,ZS_updated,XS_updated,REFC_updated,REFS_updated,RigidLift_updated,ElasticLift_updated,ElasticMoment_updated] = remove_duplicated_YS(YS,ZS,XS,REFC,REFS,RigidLift,ElasticLift,ElasticMoment);

