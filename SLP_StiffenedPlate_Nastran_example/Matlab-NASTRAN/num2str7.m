function str = num2str7(num)


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


str(8)=' ';

end
