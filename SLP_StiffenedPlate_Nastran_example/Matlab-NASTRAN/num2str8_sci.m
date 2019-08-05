function str = num2str8_sci(num)


if num<1e-4 && num>0
    
    str = '        ';
    
    str_temp = num2str(num);
    
    length_str_temp = length(str_temp);
    
    if length_str_temp>8
        
        last_4_digit = str_temp(end-3:end);
        first_4_digit = str_temp(1:4);
        
        str=[first_4_digit last_4_digit];
        
    else
        
        str(1:length_str_temp) = str_temp;
        
    end
    
else
    
    
    str = num2str8(num);
end