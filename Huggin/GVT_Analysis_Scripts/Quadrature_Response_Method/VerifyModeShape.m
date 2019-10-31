function VerifyModeShape(Info,FreqDomain,ModeFreqHz,ModeShape,fh)
% function VerifyModeShape(FreqDomain,ModeFreqHz,ModeShape,fh)
%
% This function generates frequency response plots of the raw data
% and then overlays the magnitude/phase of the estimated mode shapes.
%
% Inputs:
%   Info, FreqDomain - Data structures returned by GVT_Load.
%   ModeFreqHz - Nm-by-1 vector of estimated modal frequencies in Hz.
%   ModeShape - Naccel-by-Nm vector of mode shapes. The i^th column,
%             ModeShape(:,i), contains the i^th mode shape.
%   fh - Optional input to specify the plot figure handle [Default: fh=1].
%
% Outputs:
%   None

% Input parsing
nin=nargin;
if nin==4
    fh = 1;
end

% Grab experimental frequency response data
Accel2Run = Info.Accel2Run;
G = FreqDomain.FreqResponse;
fHz = FreqDomain.FrequencyHz;
Naccel = size(G(:,:),2);
Nmode = numel(ModeFreqHz);

% Zoom to frequency of interest
FreqWinHz = [ModeFreqHz(1)-1 ModeFreqHz(end)+1];

% Smoothed freq response in window of interest
% This is used to generate a smoother phase plot.  The magnitude
% plot is not smoothed as the peaks are critical in the mode id.
Gs = G;
Nw = 5;
fidx = find( fHz>FreqWinHz(1) & fHz<FreqWinHz(end));
Gs(fidx,:) = mysmooth(G(fidx,:),Nw);

% Plot frequency response results for each accel
for j=1:Naccel
    % Estimated mode shape gain and phase
    ModeShapeGain = 20*log10(abs(ModeShape(j,:)));
    ModeShapePhase = 90*ones(Nmode,1);
    ModeShapePhase( ModeShape(j,:) <0 ) = -90;
    
    % Get mag/phase corresponding to accel j
    run = Accel2Run(j,1);
    meas = Accel2Run(j,2);
    Gmag = 20*log10(abs(G(:,run,meas)));
    Gph = angle(G(:,run,meas))*180/pi;
    
    % Smoothed phase and avoid wrapping at 180degs
%     Gphs = angle(Gs(:,run,meas));
%     for i=1:Nmode
%         fidxi = find( fHz>ModeFreqHz(i)-0.25 & fHz<ModeFreqHz(i)+0.25);
%         Gphs(fidxi) = unwrap(Gphs(fidxi));
%     end
%     Gphs = Gphs*180/pi;        
    Gphs = angle(Gs(:,run,meas))*180/pi; 
    wrapangle = 25;            % Plot from [-180 180]-wrapangle
    idx = find(Gphs>180-wrapangle); 
    Gphs(idx)=Gphs(idx)-360;
    
    % Show Results
    figure(fh);
    
    subplot(2,1,1)
    semilogx(fHz,Gmag,'b',ModeFreqHz,ModeShapeGain,'ro');
    ylabel('|G(jw)| in dB')
    title([' Accel # = ' int2str(j)]);
    xlim(FreqWinHz);
    
    idx = find(fHz>FreqWinHz(1) & fHz<FreqWinHz(2));
    ymax = max( max(ModeShapeGain), max(Gmag(idx)) );
    ymin = min( min(ModeShapeGain), min(Gmag(idx)) );
    ylim([ymin ymax]);
    
    subplot(212)
    semilogx(fHz,Gph,'b:',fHz,Gphs,'g',ModeFreqHz,ModeShapePhase,'ro');
    xlabel('Frequency (Hz)')
    ylabel('\angle G(jw) in degs')
    xlim(FreqWinHz);
    %ylim([-180 180]-wrapangle);
    
    % Print Results
    fprintf('\nAccel # = %d \n',j);
    for i=1:Nmode
       fprintf('Mode = %d, Freq = %4.3f Hz, Amp=%4.3f, \n',...
           i,ModeFreqHz(i),ModeShape(j,i))
    end
    
    
    pause
end

