function write_SPC1(design_folder,bdf_file_name, SET1_nodes,set1_label,start,cardname,align_line)
% This program is to set up the spline4 card in NASTRAN at given structural
% nodes
%

% set1_label
% give node label, write them into NASTRAN format set1

% addpath([PATH.main_program '\subroutines\nastran_str']);

% set1_label=SET1.set1_label;

% when using SPC1, one needs to add component (1-6) to the node list
% ==============================
set1_nodes_number=length(SET1_nodes); % number of nodes in SET1

% set1_bdf_lines_num=floor((set1_nodes_number-7)/8)+1; % number of bdf lines in SET1


if set1_nodes_number<=7
    set1_bdf_lines_num=1;
    
elseif set1_nodes_number>7
    
    set1_bdf_lines_num=ceil((set1_nodes_number-7)/8)+1;
    
end


set1_bdf_each_line=['                                                                        '];


if start==1
    fid102=fopen([design_folder filesep bdf_file_name],'w+');fclose(fid102);
    fid102=fopen([design_folder filesep bdf_file_name],'a');
else
    fid102=fopen([design_folder filesep bdf_file_name],'a');
end

fprintf(fid102,'%s\n',align_line);

%%%%% only one line
if set1_bdf_lines_num<=1
    %     line_label=1;
    % the first line
    set1_bdf_each_line(1:length(cardname))=cardname;% 1-8
    set1_bdf_each_line(9:16)=num2str_integer(set1_label);% 9-16
    
    for nodenum=1:set1_nodes_number
        
        set1_bdf_each_line((nodenum+1)*8+1:(nodenum+2)*8)=num2str_integer(SET1_nodes(nodenum));
        
    end
    
    
    fprintf(fid102,'%s\n',set1_bdf_each_line);
    
    %%% more than one line
    
elseif set1_bdf_lines_num>1
    
    line_label=1;
    % the first line
    set1_bdf_each_line(1:length(cardname))=cardname;% 1-8
    set1_bdf_each_line(9:16)=num2str_integer(set1_label);% 9-16
    
    for nodenum=1:7
        
        set1_bdf_each_line((nodenum+1)*8+1:(nodenum+2)*8)=num2str_integer(SET1_nodes(nodenum));
        
    end
    
    fprintf(fid102,'%s\n',set1_bdf_each_line);
    
    %     line_set1_bdf{set1_bdf_line_num+line_label}=set1_bdf_each_line;
    
    residual_node_num=(set1_nodes_number-7)-(1-1)*8;
    
    line_label=line_label+1;
    
    while residual_node_num>=8
        
        set1_bdf_each_line(1:8)=['        '];
        
        for node_num_line=1:8
            set1_bdf_each_line(node_num_line*8+1:node_num_line*8+8)=...
                num2str_integer(SET1_nodes(7+8*(line_label-2)+node_num_line));
            
        end
        fprintf(fid102,'%s\n',set1_bdf_each_line);
        %         line_set1_bdf{set1_bdf_line_num+line_label}=  set1_bdf_each_line;
        residual_node_num=(set1_nodes_number-7)-(line_label-1)*8;
        
        line_label=line_label+1;
        
    end
    
    set1_bdf_each_line=['                                                                        '];
    for node_num_line=1:residual_node_num
        set1_bdf_each_line(node_num_line*8+1:node_num_line*8+8)=...
            num2str_integer(SET1_nodes(7+8*(line_label-2)+node_num_line));
    end
    
    %     line_set1_bdf{set1_bdf_line_num+line_label}=  set1_bdf_each_line;
    fprintf(fid102,'%s\n',set1_bdf_each_line);
end


fclose(fid102);

