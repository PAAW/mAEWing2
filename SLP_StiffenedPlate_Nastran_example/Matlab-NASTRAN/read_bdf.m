function [bdf_elem,bdf_node,bdf_cord]=read_bdf(filename)

% this function is to read nodes and element connectivity and save in one
% file bdf_node_elem
% Revise the read of scientific data; March-19-2016
bdf_elem=[];
bdf_node=[];
bdf_cord=[];
fid100=fopen(filename);

if fid100>0
    % string1='CQUAD4';string2='GRID';string3='CBAR';
    
    status=fseek(fid100,0,'eof');
    EOF=ftell(fid100);
    currentFPI=fseek(fid100,0,'bof');
    
    
    elem_num=1;
    node_num=1;
    cord_num=1;
    
    bdf_cord=zeros(1,10);
    
    while currentFPI<EOF
        
        linebdf=fgetl(fid100);currentFPI=ftell(fid100);
        
        if length(linebdf)>=7
            %
            if strcmp(linebdf(1:6),'CQUAD4')
                
                bdf_elem_bdf{elem_num}=linebdf;
                bdf_elem_temp=linebdf;
                bdf_elem(elem_num,:)=str2num(bdf_elem_temp(9:end));
                
                elem_num=elem_num+1;
                
            elseif strcmp(linebdf(1:6),'CTRIA3')
                
                bdf_elem_bdf{elem_num}=linebdf;
                bdf_elem_temp=linebdf;
                bdf_elem(elem_num,:)=str2num(bdf_elem_temp(9:end));
                
                elem_num=elem_num+1;
            end
            
            
            %
            if strcmp(linebdf(1:4),'GRID')
                bdf_node(node_num,1)=str2num(linebdf(9:16)); % node ID
                
                
                %% X
                 if strfind(linebdf(25:32),'-')% THE FIRST negative sign could be a negative value
                    
                    %                 disp('============= yes =============');
                    negative_loc=strfind(linebdf(25:32),'-');
                    
                    bdf_node(node_num,2)=str2num(linebdf(25:24+ negative_loc-1))*10^(-str2num(linebdf(24+ negative_loc+1:32)));
                    
                else
                    
                    bdf_node(node_num,2)=str2num(linebdf(25:32)); % node Y cord
                    
                 end
                 
                 
                %% Y
                if strfind(linebdf(33:40),'-')% THE FIRST negative sign could be a negative value
                    
                    %                 disp('============= yes =============');
                    negative_loc=strfind(linebdf(33:40),'-');
                    
                    bdf_node(node_num,3)=str2num(linebdf(32+1:32+ negative_loc-1))*10^(-str2num(linebdf(32+ negative_loc+1:40)));
                    
                    
                else
                    
                    bdf_node(node_num,3)=str2num(linebdf(33:40)); % node Y cord
                    
                end
                
                
                bdf_node(node_num,4)=str2num(linebdf(41:end)); % node Z cord
                node_num=node_num+1;
            end
            %
            
            % read coordinate frames 302 (302 is updated as the spars and ribs locations change)
            if strcmp(linebdf(1:6),'CORD2R') && str2double(linebdf(9:16))==302
                
                bdf_cord(cord_num,1)=str2double(linebdf(9:16));
                
                bdf_cord(cord_num,2)=str2double(linebdf(25:32));
                bdf_cord(cord_num,3)=str2double(linebdf(33:40));
                bdf_cord(cord_num,4)=str2double(linebdf(41:48));
                
                bdf_cord(cord_num,5)=str2double(linebdf(49:56));
                bdf_cord(cord_num,6)=str2double(linebdf(57:64));
                bdf_cord(cord_num,7)=str2double(linebdf(65:end));
                
                % next line
                linebdf=fgetl(fid100);currentFPI=ftell(fid100);
                bdf_cord(cord_num,8)=str2double(linebdf(9:16));
                bdf_cord(cord_num,9)=str2double(linebdf(17:24));
                bdf_cord(cord_num,10)=str2double(linebdf(25:end));
                
                cord_num=cord_num+1;
                
            end
        end
        
    end
    
else
    bdf_elem=NaN;
    bdf_node=NaN;
    disp('2D planform mesh failed!!!!!!!!!!!!!!!!!!')
end
