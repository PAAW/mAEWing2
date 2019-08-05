function trim_var_value = read_trim_var(sol144_f06,trim_variable_name,cases_num)

trim_variable_name{1}='INTERCEPT';

% trim_variable_name={
%     'INTERCEPT'
%     'OWLEADL'
%     'IWBODYL'
%     'OWLEADR'
%     'IWBODYR'
%     'ANGLEA'
%     'PITCH'
%     'URDD3'
%     'URDD5'
%     'OWTEL1'
%     'OWTEL2'
%     'OWTEL3'
%     'OWTEL4'
%     'OWTER4'
%     'OWTER3'
%     'OWTER2'
%     'OWTER1' };

trim_variables_number =length(trim_variable_name);


trim_var_value = zeros( trim_variables_number,1);


fid=fopen(sol144_f06);

status=fseek(fid,0,'eof');
EOF=ftell(fid);
currentFPI=fseek(fid,0,'bof');


linedata=fgetl(fid);currentFPI=ftell(fid); % first line of sol144.f06




keyword = 'AEROELASTIC TRIM VARIABLES';
% TRIM_Flag = 1;

case_num = 0;

while currentFPI<EOF && case_num<=cases_num
    
    
    linedata=fgetl(fid);currentFPI=ftell(fid);
    
    str = findstr(linedata,keyword);
    
    if length(linedata)==80 && isempty(str) == 0
        case_num = case_num+1;
        trim_var_num = 1;
        linedata=fgetl(fid);currentFPI=ftell(fid);
        linedata=fgetl(fid);currentFPI=ftell(fid);
        linedata=fgetl(fid);currentFPI=ftell(fid);
        
        % fourth -line after 'AEROELASTIC TRIM VARIABLES';
        linedata=fgetl(fid);currentFPI=ftell(fid);
        
        %         indexC = strfind(linedata,trim_variable_name(trim_var_num));
        TRIM_VAR = trim_variable_name{trim_var_num};
        
        TRIM_exist = findstr(linedata,TRIM_VAR);
        
        if isempty(TRIM_exist)==0
            
            TRIM_Flag = 1;
        else
            TRIM_Flag = 0;
        end
        
        
        while  TRIM_Flag && trim_var_num<=trim_variables_number% && isempty(indexC) == 0
            
            trim_var_value(trim_var_num,case_num) = str2num(linedata(92:104));
            
            
            trim_var_num=trim_var_num+1;
            
            
            if TRIM_Flag  && trim_var_num <= trim_variables_number
                
                linedata=fgetl(fid);currentFPI=ftell(fid);
                
                TRIM_VAR = trim_variable_name{trim_var_num};
                
                TRIM_exist = findstr(linedata,TRIM_VAR);
                
                if isempty(TRIM_exist)==0
                    
                    TRIM_Flag = 1;
                else
                    TRIM_Flag = 0;
                end
                
                %                 indexC = strfind(linedata,trim_variable_name(trim_var_num));
            end
            
        end
        
    end
    
end

fclose('all');
