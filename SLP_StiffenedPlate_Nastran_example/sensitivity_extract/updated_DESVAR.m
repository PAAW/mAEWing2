function updated_DESVAR(x0)


filename = [pwd filesep 'verification_stiffened_plate\DESVAR.dat'];

fid111=fopen(filename);
status=fseek(fid111,0,'eof');
EOF=ftell(fid111);
currentFPI=fseek(fid111,0,'bof');

line_num=1; % line_num should be consistent with DESVAR number

while currentFPI<EOF
    
    linef06=fgetl(fid111);currentFPI=ftell(fid111);
    
    desvar_dat{line_num}=linef06;
    
    desvar_temp = desvar_dat{line_num};
    
    desvar_temp(25:32) = num2str8(x0(line_num));
    desvar_dat{line_num} = desvar_temp;
    line_num=line_num+1;
    
    
end


%% write files
new_filename = [pwd filesep 'verification_stiffened_plate\DESVAR_updated.dat'];
fid102=fopen(new_filename,'wt');

for jjj=1:line_num-1
    
    fprintf(fid102,'%s\n',[desvar_dat{jjj}]);
    
end


fclose('all');




