function string=nastran_str_number(number,str)

% combine str and number for NASTRAN string


length_of_num = length(num2str(number));

length_of_str = length(str);

string='        ';


if length_of_num+length_of_str>8
    
    disp('---reduce string ---');
    
    
    
else
    
    
    string(1:length_of_str)=str;
    
    string(length_of_str+1:length_of_num+length_of_str) = num2str(number);
    
end




