function write_FORCE(nodelabel,forceID,force,design_folder,fname,vector,start)


if strcmp(start,'1')
fid102=fopen([design_folder filesep fname],'w+');fclose(fid102);
fid102=fopen([design_folder filesep fname],'a');

else
    
fid102=fopen([design_folder filesep fname],'a');
end

for ii = 1:length(nodelabel)
    
    force_sample = ['FORCE    3       21      0       500.   ' num2str8(vector(1)) num2str8(vector(2)) num2str8(vector(3))];
    
    force_sample(9:16) = num2str_integer(forceID);
    
    force_sample(17:24) = num2str_integer(nodelabel(ii));
    
    force_sample(33:40) = num2str8(force(ii));
    
    fprintf(fid102,'%s\n',force_sample);
    
end

fclose('all');



