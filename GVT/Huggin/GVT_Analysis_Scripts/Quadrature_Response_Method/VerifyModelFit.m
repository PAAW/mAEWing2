function VerifyModelFit(Info,FreqDomain,fHzfit,Gfit,fh)
% function VerifyModelFit(Info,FreqDomain,fHzfit,Gfit,fh)
%
% This function generates frequency response plots of the raw data
% and then overlays the magnitude/phase of the estimated mode shapes.
%
% Inputs:
%   Info, FreqDomain - Data structures returned by GVT_Load.
%   FitFreqHz - Nw-by-1 vector of frequencies in Hz.
%   Gfit - Frequency response fit. This should have the same dimension as 
%     FreqDomain.FreqResponse.
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

% Zoom to frequency of interest
FreqWinHz = [fHzfit(1) fHzfit(end)];

% Smoothed freq response in window of interest
% This is used to generate a smoother phase plot.  The magnitude
% plot is not smoothed as the peaks are critical in the mode id.
Gs = G;
Nw = 5;
fidx = find( fHz>FreqWinHz(1) & fHz<FreqWinHz(end));
Gs(fidx,:) = mysmooth(G(fidx,:),Nw);

% Plot frequency response results for each accel
for j=1:Naccel    
    % Get mag/phase corresponding to measured accel j
    run = Accel2Run(j,1);
    meas = Accel2Run(j,2);
    Gmag = 20*log10(abs(G(:,run,meas)));
    Gph = angle(G(:,run,meas))*180/pi;
    
    % Get mag/phase corresponding to model fit for accel j
    Gfitmag = 20*log10(abs(Gfit(:,run,meas)));
    Gfitph = angle(Gfit(:,run,meas))*180/pi;
    
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
    semilogx(fHz,Gmag,'b',fHzfit,Gfitmag,'r.');
    ylabel('|G(jw)| in dB')
    title([' Accel # = ' int2str(j)]);
    xlim(FreqWinHz);
    
    idx = find(fHz>FreqWinHz(1) & fHz<FreqWinHz(2));
    ymax = max( max(Gfitmag), max(Gmag(idx)) );
    ymin = min( min(Gfitmag), min(Gmag(idx)) );
    ylim([ymin ymax]);
    
    subplot(212)
    semilogx(fHz,Gph,'b:',fHz,Gphs,'g',fHzfit,Gfitph,'r.');
    xlabel('Frequency (Hz)')
    ylabel('\angle G(jw) in degs')
    xlim(FreqWinHz);
    %ylim([-180 180]-wrapangle);
    
    
    pause
end

