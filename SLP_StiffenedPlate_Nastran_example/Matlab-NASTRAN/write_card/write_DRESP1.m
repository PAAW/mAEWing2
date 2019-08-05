function write_DRESP1(ID,response_label,RTYPE,component,element_node_id,...
    design_folder,file_name,start)
% not finished code


line_sample ='DRESP1  11      D3      DISP                    3               199';

if strcmp(start,'1')
    fid101=fopen([design_folder filesep file_name],'w+');fclose(fid101);
    fid101=fopen([design_folder filesep file_name],'a');
else
    fid101=fopen([design_folder filesep pfile_name],'a');
end


