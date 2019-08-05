function write_PCOMP_v1(design_folder,pcomp_bdf_file_name,...
                               PCOMP_label,section,...
                               MatID,Layerthickness,FiberPlyAngle,...
                                Sym_flag,Z0,start,FT)


% write the file with given layer thickness, offset, sym or not.
% PCOMP  -- PID -- Z0 -- NSM -- SB -- FT -- TREF -- GE -- LAM
% Z0 - Distance from the reference plane to the bottom surface
%
%
if strcmp(start,'1')
    fid101=fopen([design_folder filesep pcomp_bdf_file_name],'w+');fclose(fid101);
    fid101=fopen([design_folder filesep pcomp_bdf_file_name],'a');
else
    fid101=fopen([design_folder filesep pcomp_bdf_file_name],'a');
end


pcomp_line_dat = ['$ ####  ======= PCOMP for ' section '  ============== ##########'];
fprintf(fid101,'%s\n',pcomp_line_dat);


if isempty(Z0)
    Z0 = blanks(8);
else
    Z0 = num2str8(Z0);
    
end

if strcmp(Sym_flag,'SYM')  
    pcomp_line_dat = ['PCOMP,' num2str_integer(PCOMP_label) ','  Z0 ', ,'  FT ' ,,               SYM'];
    fprintf(fid101,'%s\n',pcomp_line_dat);

else

    pcomp_line_dat = ['PCOMP,   ' num2str_integer(PCOMP_label) ','   Z0 ',   ,             '  FT];
    fprintf(fid101,'%s\n',pcomp_line_dat);
 
end



%% 

for layer_num_temp = 1:length(Layerthickness)
    
        
    pcomp_line_dat =['        ' num2str_integer(MatID(layer_num_temp ))  num2str8(Layerthickness(layer_num_temp))  num2str8(FiberPlyAngle(layer_num_temp)) 'YES'];
    fprintf(fid101,'%s\n',pcomp_line_dat);
    

end


%%
fclose(fid101);
