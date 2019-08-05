function string=nastran_str_number_trim(number,str)

% combine str and number for NASTRAN string
% not asked to be 8-bit string as output
% referred in 
%             trim_variable_name = NASTRAN_Trim_Var(TE_flaps_num)

length_of_num = length(num2str(number));

length_of_str = length(str);

% string='        ';


if length_of_num+length_of_str>8
    
    disp('---reduce string ---');
    
    
    
else
    
    
    string(1:length_of_str)=str;
    
    string(length_of_str+1:length_of_num+length_of_str) = num2str(number);
    
end




