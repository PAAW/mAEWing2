function write_CQUAD4(elements,design_folder,fname)



fid102=fopen([design_folder filesep fname],'w+');fclose(fid102);
fid102=fopen([design_folder filesep fname],'a');



for ii = 1:size(elements,1)
    
    cquad4 ='CQUAD4   1       1      1       2       23      22      0';
    
    cquad4(9:16) = num2str_integer(elements(ii,1));
    cquad4(17:24) = num2str_integer(elements(ii,2));
    
    cquad4(25:32) = num2str_integer(elements(ii,3));
    cquad4(33:40) = num2str_integer(elements(ii,4));
    cquad4(41:48) = num2str_integer(elements(ii,5));
    cquad4(49:56) = num2str_integer(elements(ii,6));
    
    
    fprintf(fid102,'%s\n',cquad4);
    
end
fclose('all');
