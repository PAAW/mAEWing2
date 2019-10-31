function [] = GVT_PlotModeShapes(Info,Modal,ModeShapes)

Nmode = Modal.Nmode;
Nrun = Info.Nrun;
% ForwardAccel = Info.ForwardAccel;
% AftAccel = Info.AftAccel;
AccelPos = Info.AccelPos;
xcoord = AccelPos(:,1);             % X-coordinates of accelerometers
ycoord = AccelPos(:,2);             % Y-coordiantes of accelerometers

PlotOrder = [11:-2:1 29 25 27 13:2:23 24:-2:14 28 26 30 2:2:12];

%% Plot mode shapes

% 3-D plots
for i = 1:Nmode
    zcoord=ModeShapes(:,i);
    
    figure
    title(['Exact mode shape for mode:',num2str(i),...
        ', Freq:',num2str(Modal.w_mode(i)),'rad/s']);
    hold on;
    p=fill3(xcoord(PlotOrder),ycoord(PlotOrder),zcoord(PlotOrder),'b');
    set(p,'EdgeColor','b')
    view(10,20);
    daspect([max(daspect)*[1 1] 1])
end


% 2-D plots
% heave = zeros(Nmode,Nrun);
% twist = zeros(Nmode,Nrun);
% for i = 1:Nmode
%     
%     heave(i,:) = (2.75 *ModeShapes(ForwardAccel,i) + ...
%                 3.25 * ModeShapes(AftAccel,i)) ... 
%                 /(2*(2.75 + 3.25));
%     twist(i,:) = atand( (ModeShapes(ForwardAccel,i) - ...
%                 ModeShapes(AftAccel,i))...
%                     /6);
%     figure
%     hold on
%     plot(AccelPos(ForwardAccel,1),heave(i,:)*4,'b')
%     plot(AccelPos(ForwardAccel,1),twist(i,:),'r')
%     legend('Heave','Twist')
%     title(['Mode shape for mode:',num2str(i),...
%         ', Freq:',num2str(Modal.w_mode(i)),'rad/s']);
%     grid on
%     garyfyFigure
end