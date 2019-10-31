%% Read_GVT_Data.m: Uploads and reads the GVT data

%% Clear workspace
clear
close all
clc


%% Load GVT Data
load GVT_Data_Huggin.mat


%% GVT Parameters

% Position of accelerometers
% AccelPos = [x,y]
% x = positive towards left wing, y = positive towards aft
AccelPos = Info.AccelPos; 

% Position of the shaker
ShakerPos = Info.ShakerPos;

% Frequency Response
G = FreqDomain.FreqResponse;        % Frequency response data from GVT
fHz = FreqDomain.FrequencyHz;       % Corresponding frequencies in Hz

FreqWin_Hz = Info.FreqWin_Hz;       % The frequency range of GVT


%% Accelerometer Positions
open Accel_Pos.fig

%% Obtain frequency response for a particular accelerometer

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose the accelerometer in the line below
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
accel_num = 10;


% Obtain the index of the frequency response data
Run2Accel = Info.Run2Accel;            
[idx1,idx2]=find(Run2Accel==accel_num);  

% Frequency response of the particular accelerometer
G_accel_num = G(:,idx1,idx2);

% Plot the magnitude of the data
figure
Gmag = 20*log10(abs(G_accel_num));
semilogx(fHz,Gmag,'b');
xlabel('Frequency (Hz)')
ylabel('|G(jw)| in dB')
title(['Accel # = ' int2str(Run2Accel(idx1,idx2))]);
xlim(FreqWin_Hz);
ylim([-80 20]);
grid minor