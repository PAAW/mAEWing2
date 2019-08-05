test_case_element_connectivity =   Elements_Label(:,2:end);

test_case_node_cords = GRID_cord;

for elem = 1:size(test_case_element_connectivity)
    
    
    
    node_labels = test_case_element_connectivity(elem,:);
    
    
    % center node coordinate
    
    
    centerX = sum( test_case_node_cords(node_labels,2))/length(node_labels);
    centerY = sum( test_case_node_cords(node_labels,3))/length(node_labels);
    
    for layer = 1:size(T0T1,1)
        
        T0 = T0T1(layer,1);
        T1 = T0T1(layer,2);
        
        theta(layer) = VAT_fiber_ply_angle_1D(T0,T1,centerX,center,width);
        
        
    end
    
    
% % % %     figure(201);hold on;plot(centerX,centerY,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0],'Marker','o',...
% % % %     'LineStyle','none',...
% % % %     'Color',[0 0 0])
% % % %     
% % % %     if abs(theta(layer))<1e-5
% % % %       u = 1; v=0;
% % % %     else
% % % %        u= cosd(theta(layer));
% % % %        v =sind(theta(layer));
% % % %     end
% % % %     quiver(centerX,centerY,u,v,'AutoScaleFactor',0.05,'Color',[0 0 0])
% % % %     
% % % %     
    %% write PCOMP
    
    section  = ['VAT=== Elem#' num2str(elem)];
    Sym_flag = 'ASYM';
    FT = ' ';
    PCOMP_label = elem;
    
    MatID          = 777;
    Layerthickness = 1e-3*ones(length(theta),1);
    Z0 ='';
    FullPlyAngles = theta;
    
    if elem == 1
        start ='1';
    else
        start ='0';
    end
    design_folder ='C:\Users\weizhao\Documents\Paper\AIAAJ_curve_fiber\AIAAJ_Vibration_Code_modified\NASTRAN_examples\design_folder';
    pcomp_bdf_file_name = 'test_VAT_para.dat';
    write_PCOMP_v1(design_folder,pcomp_bdf_file_name,...
        PCOMP_label,section,...
        MatID,Layerthickness,FullPlyAngles,...
        Sym_flag,'',start,FT);
    
    
end