function Info = GVT_Info()
%% Outputs GVT parameters
% This function contains the parameters of the experiment and outputs them
% as a structure Info
% 
% Info has following fields
%   - Fs:= Sampling frequency of the experiment
%   - Run2Accel:= run-by-2 matrix with Run2Accel(i,j) providing
%         the accelerometer number for measurement j of run i.
%   - Accel2Run: Naccel-by-2 matrix with Accel2Run(i,1:2) providing
%         the [Run #, Measurement #] for accelometer i.
%   - AccelPos:= Naccel-by-2 matrix where  AccelPos(i,1) and
%         AccelPos(i,2) are the (X,Y) positions of accel i measured
%         in inches.  The coordinates are located at the aircraft nose
%         with Y pointing aft and X pointing left (Hence Z is up).
%   - ShakerPos:= 1-by-2 matrix where ShakerPos(1) and ShakerPos(2) are
%       the (X,Y) positions of the shaker attachment point.  The 
%       coordinates are located at the aircraft nose with Y pointing aft
%       and X pointing left (Hence Z is up).
%   - Nrun := Number of runs in the experiment
%   - Nmeas := Number of measurements in each run
%   - Naccel := Number of total acceelerations measurements
%   - Force_V2lb := Conversion factor from measured voltage to force (Lbs)
%   - Accel_V2g := Conversion factor from measured voltage to acceleration
%       (g)
%   - FreqWin+Hz := Target frequency window of the experiment
%   - PlotOrder := Order of the measurement point used to plot the mode
%       shapes
%% Sample Info
Fs = 2000;                      % Sampling freq, Hz

%% Calibrations
% Force: 0.4684 V/lb , Accels: 1.028 V/g
Force_V2lb = 0.4684;
Accel_V2g = 1.028;

in2m = 0.0254;

%% Accel Mapping
% The data consists of Nrun experiments with each experiment providing
% Nmeas accelerometer measurements.  The total number of accelerometer
% measruements is thus Naccel=Nrun*Nmeas.   The accels are numbered
% from 1 to Naccel and the position of each accelerometer is provided
% in the array AccelPos below.  The mapping between accel numbering
% and experiment is provided by two variables:
%   A) Run2Accel := Nrun-by-2 matrix with Run2Accel(i,j) providing
%         the accelerometer number for measurement j of run i.
%   B) Accel2Run := Naccel-by-2 matrix with Accel2Run(i,1:2) providing
%         the [Run #, Measurement #] for accelometer i.
Run2Accel = [1 2; 3 4; 5 6; 7 8; 9 10; 11 12; 13 14; 15 16; 17 18; ...
    19 20; 21 22; 23 24; 25 26; 27 28; 29 30];
[Nrun,Nmeas] = size(Run2Accel);
Naccel= Nrun*Nmeas;
Accel2Run = zeros(Naccel,2);
for i=1:Naccel;
    [ridx,cidx]=find(Run2Accel==i);
    Accel2Run(i,:) = [ridx cidx];
end


%% Accel Positions
% As noted above, the accelerometer measurements recorded in all
% experimental runs are numbered from 1 to Naccel.  The positions
% of these accelerometers is provided by:
%     AccelPos := Naccel-by-2 matrix where  AccelPos(i,1) and
%         AccelPos(i,2) are the (X,Y) positions of accel i measured
%         in inches.  The coordinates are located at the aircraft nose
%         with Y pointing aft and X pointing left (Hence Z is up).

T = [cosd(-22) sind(-22);
    -sind(-22) cosd(-22)] ;
AccelPos(:,1:12) = repmat(in2m*[9.5 22.286]',1,12) + ...
                    T * [0.2036	-0.0619;
                        0.2543	0.1265;
                        0.5084	-0.0619;
                        0.5084	0.0651;
                        0.8957	-0.0619;
                        0.8957	0.0651;
                        1.2767	-0.0619;
                        1.2767	0.0651;
                        1.6577	-0.0619;
                        1.6577	0.0651;
                        1.9562	-0.0619;
                        1.9562	0.0651;
                        ]';
                
      
T = [cosd(22) sind(22);
    -sind(22) cosd(22)] ;
AccelPos(:,13:24) = repmat(in2m*[-9.5 22.286]',1,12) + ...
                     T*[-0.2036	-0.0619;
                        -0.2543	0.1265;
                        -0.5084	-0.0619;
                        -0.5084	0.0651;
                        -0.8957	-0.0619;
                        -0.8957	0.0651;
                        -1.2767	-0.0619;
                        -1.2767	0.0651;
                        -1.6577	-0.0619;
                        -1.6577	0.0651;
                        -1.9562	-0.0619;
                        -1.9562	0.0651;
                        ]';
                
                
AccelPos(:,25:26) = [0 13;
                     0 38.523]'*in2m;
AccelPos(:,27:28) = [-9.5 19.275;
                     -9.5 38.523]'*in2m;
AccelPos(:,29:30) = [9.5 19.275;
                     9.5 38.523]'*in2m;
AccelPos = AccelPos';

%% Shaker Position
% The shaker position remains same for all the runs for this particular
% experiment.
%       ShakerPos := 1-by-2 matrix where ShakerPos(1) and ShakerPos(2) are
%       the (X,Y) positions of the shaker attachment point.  The 
%       coordinates are located at the aircraft nose with Y pointing aft
%       and X pointing left (Hence Z is up).
ShakerPos = [52.15 35.39] * in2m;

%% Frequency Range
% The target frequency range
FreqWin_Hz = [3 35];

%% Plot order
% Order of points in the order they are to be plotted

% PlotOrder = [20:-2:12 2:2:10 9:-2:1 11:2:19];
PlotOrder = [11:-2:1 29 25 27 13:2:23 24:-2:14 28 26 30 2:2:12];
%% Store Experiment Info
Info.Fs = Fs;
Info.Run2Accel = Run2Accel;
Info.Accel2Run = Accel2Run;
Info.AccelPos = AccelPos;
Info.ShakerPos = ShakerPos;
Info.Nrun = Nrun;
Info.Nmeas = Nmeas;
Info.Naccel = Naccel;
Info.Force_V2lb = Force_V2lb;
Info.Accel_V2g = Accel_V2g;
Info.FreqWin_Hz = FreqWin_Hz;
Info.PlotOrder = PlotOrder;
end