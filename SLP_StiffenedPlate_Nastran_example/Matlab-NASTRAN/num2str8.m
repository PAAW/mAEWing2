function str = num2str8(num)


% GENERATE STR FOR NASTRAN INPUT FILES

string =  deblank(num2str(num,'%g'));


ind_e = find(string=='e');
ind_d = find(string=='.');

if (isempty(ind_d) && length(string)<8)
    if (isempty(ind_e))
        string = [string,'.'];
    else
        string = [string(1:ind_e-1),'.',string(ind_e:end)];
    end
    str = [string,  blanks(8-length(string))];
    return
end
%
if (isempty(ind_e))
    nm = min(8, length(string));
    str = [string(1:nm), blanks(8-nm)];
else
    le = length(string) +1 - ind_e;
    str = [string(1:(length(string)-le)), string(ind_e:end)];
end

str = [str,  blanks(8-length(str))];

if length(str)>8
    
    str=sprintf('%8f',num);
    
    if length(str)>8
        
        str=str(1:8);
    end
    
end


%% small value

if num<=1e-3 && num>0
    
    
    
    string =  deblank(num2str(num,'%g'));
    
    if length(string)>8
        
        
        ind_e = find(string=='e');
        ind_d = find(string=='.');
        
        ind_zero = find(string=='0');
        
        
        num_string = string(ind_zero(end)+1:end);
        
        sci_string = ['-' num2str(ind_zero(end)-1)];
        
        
        left_num = 8-length(sci_string);
        
        
        num_string_updated = [num_string(1) '.' num_string(2:end) ];
        
        if length(num_string_updated) + length(sci_string)>8
            
            str =[num_string_updated(1: left_num) sci_string];
            
            
        else
            
            str =[num_string_updated sci_string blanks(8-length(num_string_updated)-length(sci_string))];
            
            
            
        end
    end
end

