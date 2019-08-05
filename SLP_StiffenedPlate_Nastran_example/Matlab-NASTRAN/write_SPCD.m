function write_SPCD(nodelabel,SPCD_id,deltaU,design_folder,fname,component,start)
% when use this for 

if strcmp(start,'1')
    fid102=fopen([design_folder filesep fname],'w+');fclose(fid102);
    fid102=fopen([design_folder filesep fname],'a');
    
else
    
    fid102=fopen([design_folder filesep fname],'a');
end

for ii = 1:length(nodelabel)
    
    
    spcd_sample ='SPCD     1       121     1      -1.-4   ';
    
    %     spcd_sample = ['SPCD     3       21      0       500.   ' num2str8(component(1)) num2str8(component(2)) num2str8(component(3))];
    
    
    for comp = 1:length(component)
        
        spcd_sample(9:16) = num2str_integer(SPCD_id);
        spcd_sample(17:24) = num2str_integer(nodelabel(ii));
        
        spcd_sample(25:32) = num2str_integer(component(comp));
        
        
        spcd_sample(33:40) = num2str8(deltaU);
        
        fprintf(fid102,'%s\n',spcd_sample);
        
    end
    
end

fclose('all');