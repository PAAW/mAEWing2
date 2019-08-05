function generate_RBE2_xAEWing3(design_folder, RBE2_bdf_file,...
                                dependent_label, independent_label,...
                                start, RBE2_ID,RBE2_ID_offset,node_offset,...
                                message_line,dofs)

%

if strcmp(start,'1')
    
    fid101=fopen([design_folder filesep RBE2_bdf_file],'w+');fclose(fid101);
    fid102=fopen([design_folder filesep RBE2_bdf_file],'a+');
elseif strcmp(start,'0')
    fid102=fopen([design_folder filesep RBE2_bdf_file],'a+');
end

fprintf(fid102,'%s\n',message_line);
fprintf(fid102,'%s\n','$1......2.......3.......4.......5.......6.......7.......8.......');


RBE2_line = ['RBE2    ' num2str_integer(RBE2_ID) num2str_integer(independent_label) num2str_integer(dofs) num2str_integer(dependent_label)];
fprintf(fid102,'%s\n',RBE2_line);


RBE2_line = ['RBE2    ' num2str_integer(RBE2_ID+RBE2_ID_offset) num2str_integer(independent_label+node_offset) ...
    num2str_integer(dofs) num2str_integer(dependent_label+node_offset)];
fprintf(fid102,'%s\n',RBE2_line);


fclose(fid102);