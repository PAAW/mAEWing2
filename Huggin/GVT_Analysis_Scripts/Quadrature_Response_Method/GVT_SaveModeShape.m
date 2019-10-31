function[] = GVT_SaveModeShape(Info,Modal,ModeShapes)

Nmode = Modal.Nmode;                % Number of Modes
ForwardAccel = [11:-2:1 29 25 27 13:2:23];   % Index of forward accelerometers
AftAccel = [12:-2:2 30 26 28 14:2:24];           % Index of aft accelerometers
AccelPos = Info.AccelPos;           % Position coordinates of accels
Naccel = Info.Naccel;               % Total number of acclels
xcoord = AccelPos(:,1);             % X-coordinates of accelerometers
ycoord = AccelPos(:,2);             % Y-coordiantes of accelerometers

PlotOrder = [fliplr(AftAccel) ForwardAccel];    % Plot order for modeshape

xFor = AccelPos(ForwardAccel,1);    % x-coordinates of forward accels
yFor = AccelPos(ForwardAccel,2);    % y-coordinates of forward accels
xAft = AccelPos(AftAccel,1);        % x-coordinates of aft accels
yAft = AccelPos(AftAccel,2);        % y-coordinates of aft accles



%% Save GIF

for i = 1:Nmode
    zcoord=ModeShapes(:,i);         % Current mode shapes
    
    % File name of the saved graphic
    filename=['ModeShape:' num2str(i) '.gif'];
    
%     title(['Mode shape, ',num2str(i),...
%         ', Freq:',num2str(Modal.w_mode(i)/2/pi),'Hz']);
%     
    
    % Creating different m
    for n=[1:-0.4:-1 -.8:0.2:.8]
    ph=figure(i+100);
    title(['Mode shape, ',num2str(i),...
        ', Freq:',num2str(Modal.w_mode(i)/2/pi),'Hz']);
    
    

    plot3(xcoord(PlotOrder),ycoord(PlotOrder),0*zcoord(PlotOrder),'k.','MarkerSize',18);
    hold on;
    plot3(xcoord(PlotOrder),ycoord(PlotOrder),0*zcoord(PlotOrder),'k--','LineWidth',2);
    plot3([xcoord(PlotOrder(end)) xcoord(PlotOrder(1))],...
        [ycoord(PlotOrder(end)) ycoord(PlotOrder(1))],...
        0*[zcoord(PlotOrder(end)) zcoord(PlotOrder(1))]...
        ,'k--','LineWidth',2);
    for j=1:length(xFor)
        plot3([xAft(j) xFor(j)],[yAft(j) yFor(j)],...
        [0 0],'k--','LineWidth',2);
    end
    
        
    
    plot3(xcoord(PlotOrder),ycoord(PlotOrder),n*zcoord(PlotOrder),'r.','MarkerSize',18);
    for j=1:length(xFor)
        plot3([xAft(j) xFor(j)],[yAft(j) yFor(j)],...
        [n*zcoord(AftAccel(j)) n*zcoord(ForwardAccel(j))],'b','LineWidth',2);
    end
    
    plot3(xcoord(PlotOrder),ycoord(PlotOrder),n*zcoord(PlotOrder),'b',...
        'linewidth',3);

    plot3([xcoord(PlotOrder(end)) xcoord(PlotOrder(1))],...
        [ycoord(PlotOrder(end)) ycoord(PlotOrder(1))],...
        n*[zcoord(PlotOrder(end)) zcoord(PlotOrder(1))]...
        ,'b','LineWidth',3);
    
    zlim(max(abs(zcoord))*[-2 2])
    view(15,25);
    
    drawnow
    frame=getframe(i+100);
    im=frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    
    if n == 1;
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.2);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.2);
    end
    delete(ph)
    end
end
