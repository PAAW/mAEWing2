function wirte_TEMP(nodelabel,temp,design_folder,fname,start)


if strcmp(start,'1')
fid102=fopen([design_folder filesep fname],'w+');fclose(fid102);
fid102=fopen([design_folder filesep fname],'a');

else
    
fid102=fopen([design_folder filesep fname],'a');
end

for ii = 1:length(nodelabel)
    
    force_sample = 'TEMP    1       21      0       ';
    
    force_sample(17:24) = num2str_integer(nodelabel(ii));
    force_sample(25:32) = num2str8(temp);

    
    fprintf(fid102,'%s\n',force_sample);
    
end

fclose('all');