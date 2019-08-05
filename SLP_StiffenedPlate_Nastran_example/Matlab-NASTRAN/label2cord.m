function cords = label2cord(node_label,FEM_equiv_fullnodes)

% based on nodal label to return coordinates

% cords =[label, x,y,z]




total_nodes_number = size(FEM_equiv_fullnodes,1);


try
    for num=1:total_nodes_number
        
        
        if FEM_equiv_fullnodes(num,1) == node_label
            
            
            loc = num;
            
        end
    end
    
    
    cords = FEM_equiv_fullnodes(loc,:);
    
catch
    cords =nan;
end