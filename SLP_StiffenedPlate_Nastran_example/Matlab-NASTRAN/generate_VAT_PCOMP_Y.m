% clear load

% load test_case_node_cords.txt
% load test_case_element_connectivity.txt

test_case_node_cords =   GRID_cord;

test_case_element_connectivity =  Elements(:,3:6);
% % % design_folder ='C:\Users\weizhao\Documents\Paper\AIAAJ_curve_fiber\AIAAJ_Vibration_Code_modified\NASTRAN_examples\design_folder';
    
patch_plot(test_case_element_connectivity, test_case_node_cords,201,'skin')

% % Layer - 1
% t0 = 45;
% t1 = 0;
% 
% T0T1 = [ t0 t1 ;
%     -t0 -t1;
%     t0 t1;
%     t0 t1;
%     -t0 -t1;
%     t0 t1];


center = [Stru.length Stru.width]/2;
width  = Stru.length;


% search each element, nodal coordinate to identify the center coordinate
% and then  use x to determine the fibe ply orientation
for elem = 1:size(test_case_element_connectivity)

    node_labels = test_case_element_connectivity(elem,:);

    % center node coordinate
    
    
    centerX = sum( test_case_node_cords(node_labels,2))/length(node_labels);
    centerY = sum( test_case_node_cords(node_labels,3))/length(node_labels);
    
    for layer = 1:size(T0T1,1)
        
        Theta_0 = T0T1(layer,1);
        Theta_1 = T0T1(layer,2);
        
        theta(layer) = VAT_fiber_ply_angle_1D(Theta_0,Theta_1,centerY,center(2),Stru.width);
        
        
    end
    
    
    figure(201);hold on;plot(centerX,centerY,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0],'Marker','o',...
    'LineStyle','none',...
    'Color',[0 0 0])
    
    if abs(theta(layer))<1e-5
      u = 1; v=0;
    else
       u= cosd(theta(layer));
       v =sind(theta(layer));
    end
%     quiver(centerX,centerY,u,v,'AutoScaleFactor',0.05,'Color',[0 0 0])
    
    
    %% write PCOMP
    
    section  = ['VAT=== Elem#' num2str(elem)];
    Sym_flag = 'ASYM';
    FT = ' ';
    PCOMP_label = elem;
    
    MatID          = 777;
    Layerthickness = Laminate.layer_thickness*ones(length(theta),1);
    Z0 ='';
    FullPlyAngles = theta;
    
    if elem == 1
        start ='1';
    else
        start ='0';
    end
%     design_folder ='C:\Users\weizhao\Documents\Paper\AIAAJ_curve_fiber\AIAAJ_Vibration_Code_modified\NASTRAN_examples\design_folder';
    pcomp_bdf_file_name = 'test_VAT.dat';
    write_PCOMP_v1(design_folder,pcomp_bdf_file_name,...
        PCOMP_label,section,...
        MatID,Layerthickness,FullPlyAngles,...
        Sym_flag,'',start,FT);
    
    
end



%%
% % % % % % % % % % sol103_bdf_name = 'C:\Users\weizhao\Documents\Paper\AIAAJ_curve_fiber\AIAAJ_Vibration_Code_modified\NASTRAN_examples\design_folder\sol101_VAT.bdf';
% % % % % % % % % %  nastran = 'C:\MSC.Software\MSC_Nastran\20170\bin\nast20170.exe ';
% % % % % % % % % % cd(design_folder);
% % % % % % % % % % call_NASTRAN_sol103=[nastran ' ' sol103_bdf_name];
% % % % % % % % % % system(call_NASTRAN_sol103);
% % % % % % % % % % 
% % % % % % % % % % cd(nastran_para_folder);