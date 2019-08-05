function Inplane_stress = read_stress_sol101_pch(filename,elements_num,layers_number)

% stress in punch file
% laminated stress - layer

fid100=fopen(filename);

status=fseek(fid100,0,'eof');
EOF=ftell(fid100);
currentFPI=fseek(fid100,0,'bof');


Inplane_stress = zeros(elements_num,3,layers_number);

element_num = 1;
while currentFPI<EOF
    
    linef06=fgetl(fid100);currentFPI=ftell(fid100);
    linef06=fgetl(fid100);currentFPI=ftell(fid100);
    linef06=fgetl(fid100);currentFPI=ftell(fid100);
    linef06=fgetl(fid100);currentFPI=ftell(fid100);
    
    
    linef06=fgetl(fid100);currentFPI=ftell(fid100);
    linef06=fgetl(fid100);currentFPI=ftell(fid100);
    linef06=fgetl(fid100);currentFPI=ftell(fid100);
    
    
    
    while element_num <=elements_num && currentFPI<EOF
        
        linef06=fgetl(fid100);currentFPI=ftell(fid100);
        
%         linef06
        
        element_num = str2num(linef06(1:20));
        layer_num =   str2num(linef06(21:40));
        
        stress11 =  str2num(linef06(41:55));
        stress22 =  str2num(linef06(60:73));
        
        linef06=fgetl(fid100);currentFPI=ftell(fid100);
        
        stress12 =  str2num(linef06(21:40));
        
        linef06=fgetl(fid100);currentFPI=ftell(fid100);
        linef06=fgetl(fid100);currentFPI=ftell(fid100);
        
%          [stress11 stress22 stress12]
        
        Inplane_stress(element_num,:,layer_num) = [stress11 stress22 stress12];
        
        

    end
    
    
    
end




