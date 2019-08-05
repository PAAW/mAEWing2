function nodal_displacement = read_trim_disp(sol144_f06)

%'          1000      G      3.642679E-03  -5.419594E-04   1.043052E-02
%3.742471E-04   9.782225E-04  -3.783385E-04';


nodal_displacement =[0 0 0 0 0 0 0]; % nodal ID plus 6 displacements
jj=1;
% try

fid = fopen(sol144_f06);
keyword1 = 'G';


status=fseek(fid,0,'eof');
EOF=ftell(fid);
currentFPI=fseek(fid,0,'bof');


linedata=fgetl(fid);currentFPI=ftell(fid);


while currentFPI<EOF
    
    linedata=fgetl(fid);  currentFPI=ftell(fid);
    
    %         findstr(linedata,'G')==21
    
    if length(linedata)>110 && ~isempty(findstr(linedata,'G'))
        
        if findstr(linedata,'G')==21
            
            point_id  =  str2num(linedata(1:15));
            T1 = str2num(linedata(27:41));
            T2 = str2num(linedata(42:56));
            T3 = str2num(linedata(57:71));
            R1 = str2num(linedata(72:86));
            R2 = str2num(linedata(87:101));
            R3 = str2num(linedata(102:end));
            nodal_displacement(jj,:)=[point_id,T1,T2,T3,R1,R2,R3];
            jj = jj+1;
            
        end
    end
    
    
end






% catch
%
%     nodal_displacement = [];
% end