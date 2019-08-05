function write_CTRIA3(elements,design_folder,fname)



fid102=fopen([design_folder filesep fname],'w+');fclose(fid102);
fid102=fopen([design_folder filesep fname],'a');



for ii = 1:size(elements,1)
    
    ctria3 ='CTRIA3   1       1       252     218     210      0';
    
    ctria3(9:16)  = num2str_integer(elements(ii,1));
    ctria3(17:24) = num2str_integer(elements(ii,2));
    
    ctria3(25:32) = num2str_integer(elements(ii,3));
    ctria3(33:40) = num2str_integer(elements(ii,4));
    ctria3(41:48) = num2str_integer(elements(ii,5));
%     ctria3(49:56) = num2str_integer(elements(ii,6));
    
    
    fprintf(fid102,'%s\n',ctria3);
    
end
fclose('all');
