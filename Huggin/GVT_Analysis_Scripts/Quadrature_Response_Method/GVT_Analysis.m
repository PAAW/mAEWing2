%% Processes GVT data to obtain the modal frequency and mode shapes
% This script calls the loads the GVT data using the function GVT_Load and
% processes the dat to obtain modal frequencies and mode shapes.

%% Clear workspace
clear 
close all
clc

%% Load GVT Data
dirname = '../../Data/';
[Info,TimeDomain,FreqDomain] = GVT_Load(dirname);

% Obtain data from structures
FreqWin_Hz = Info.FreqWin_Hz;       % Target frequency window for GVT
Run2Accel = Info.Run2Accel;         % Mapping from run number to accel num 
Nmeas = Info.Nmeas;                 % Number of measurement in one run
Nrun = Info.Nrun;                   % Number of runs

G = FreqDomain.FreqResponse;        % Frequency response data from GVT
fHz = FreqDomain.FrequencyHz;       % Frequencies from GVT in Hz

%% Sigma plot of all data
Gsig = sum(abs(G(:,:)).^2,2);
Gsig = 20*log10(Gsig);

% Plot the sigma plot
figure;
semilogx(fHz,Gsig,'b');
xlabel('Frequency (Hz)')
ylabel('|G(jw)| in dB')
title('Max Singular Value');
xlim(FreqWin_Hz);
ylim([-60 70]);

%% Plot individual responses (Group by related positions)
accelsets = [1 2;
             3 4;
             5 6;
             7 8;
             9 10;
             11 12;
             13 14;
             15 16;
             17 18;
             19 20;
             21 22;
             23 24;
             25 26;
             27 28;
             29 30;
             ];
[Nset,AccelPerSet] = size(accelsets);

for i=1:Nset
    for j=1:AccelPerSet
        [idx1,idx2]=find(Run2Accel==accelsets(i,j));

        figure(i+1)
        subplot(2,2,j)
        Gmag = 20*log10(abs(G(:,idx1,idx2)));
        semilogx(fHz,Gmag,'b');
        xlabel('Frequency (Hz)')
        ylabel('|G(jw)| in dB')
        title(['Accel # = ' int2str(Run2Accel(idx1,idx2))]);
        xlim(FreqWin_Hz);
        ylim([-80 20]);
        
        
        subplot(2,2,j+2)
        Gphase = angle(G(:,idx1,idx2))*180/pi;
        semilogx(fHz,Gphase,'b');
        xlabel('Frequency (Hz)')
        ylabel('Phase| in degrees')
        title(['Accel # = ' int2str(Run2Accel(idx1,idx2))]);
        xlim(FreqWin_Hz);
%         ylim([ 20]);
    end
end


%% Fit data with a Nord-order model in the target freqeuncy windoes 

f_winidx = find( fHz>FreqWin_Hz(1) & fHz < FreqWin_Hz(end) );
w_winrad = fHz(f_winidx)*2*pi;
Nord = 20;
G_fit = zeros(numel(f_winidx),Nrun,Nmeas);
for i = 1:Nrun
    for j = 1:Nmeas
        sysin = frd(G(f_winidx,i,j),w_winrad);
        sysout = fitfrd(sysin,Nord);
        G_fit(:,i,j) = freqresp(sysout,w_winrad);
    end
end

%% Frequency domain data
FreqDomain.G_fit = G_fit;
FreqDomain.w_winrad = w_winrad;

% % xxx 
% save G_fit G_fit
% % xxx

% accelsets = [1 2; 3 4; 5 6; 7 8; 9 10];
[Nset,AccelPerSet] = size(accelsets);
for i=1:Nset
    for j=1:AccelPerSet
        [idx1,idx2]=find(Run2Accel==accelsets(i,j));

        figure(i+1)
        subplot(2,2,j)
        Gmag = 20*log10(abs(G(:,idx1,idx2)));
        semilogx(fHz,Gmag,'b');
        xlabel('Frequency (Hz)')
        ylabel('|G(jw)| in dB')
        title(['Accel # = ' int2str(Run2Accel(idx1,idx2))]);
        xlim(FreqWin_Hz);
        ylim([-80 20]);
        
        
        hold on
        Gmag_fit = 20*log10(abs(G_fit(:,idx1,idx2)));
        semilogx(FreqDomain.w_winrad/2/pi,Gmag_fit,'r');
        legend('Experimental Data','Fit Model')
        grid on
        
        subplot(2,2,j+2)
        Gphase = angle(G(:,idx1,idx2))*180/pi;
        semilogx(fHz,Gphase,'b');
        xlabel('Frequency (Hz)')
        ylabel('Phase| in degrees')
        title(['Accel # = ' int2str(Run2Accel(idx1,idx2))]);
        xlim(FreqWin_Hz);
%         ylim([- 20]);

        hold on
        Gphase_fit = angle(G_fit(:,idx1,idx2))*180/pi;
        semilogx(FreqDomain.w_winrad/2/pi,Gphase_fit,'r');
        legend('Experimental Data','Fit Model')
        grid on


    end
end




%% Identified modal parameters
% Approximate modal freq from sigma plot peaks

Modal.w_mode = [5.81 8.23 12.46 15.43 19.01 24.32 26.43 33.92]*2*pi;   % For 20 order fit 
Info.FreqWin_Hz = [3 35];

Modal.Nmode = length(Modal.w_mode);

% Assume damping value
Modal.zeta = 0.03*ones(1,Modal.Nmode);        


% Obtain and plot the mode shapes
ModeShapes = GVT_ObtainModeShapes(Modal,Info,FreqDomain);
GVT_SaveModeShape(Info,Modal,ModeShapes)
save GVT_Data_Huggin FreqDomain Info Modal TimeDomain