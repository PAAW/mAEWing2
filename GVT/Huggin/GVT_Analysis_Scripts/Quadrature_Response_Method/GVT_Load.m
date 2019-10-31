function [Info,TimeDomain,FreqDomain] = GVT_Load(dirname,GVT_code)
%% Loads the GVT data & processes it to obtain frequency domain data
% This function loads the GVT data and produces frequency domain data. The 
% function first loads the experimental data by calling the function
% GVT_Info.
% 
% Inputs:
% dirname := Name of the directory where data is stored.
%
% Outputs:
% Info: A structure containing experimental information in the following
% fields
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
%   - ForwardAccel := Order of the accelerometers on the leading edge from
%      left wing tip to the right wing tip
%   - AftAccel := Order of the accelerometers on the trailing edge from the
%       left wing tip to the right wing tip
%
% TimeDomain: A structure containing time domain data in the following
% fields (Here Nw is the number of recorded data points)
%   - Time := An Nw-by-1 array of the time stamps of Nw data points
%   - InputCommand := An Nw-by-Nrun matrix where InputCommand(i,j) contains
%        the input command signal given to the shaker at data point i for
%        run j.
%   - InputForce := An Nw-by-Nrun matrix where InputForce(i,j) contains
%        the force transfered to the aircraft at data point i for run j.
%   - OutputAccel := An Nw-by-Nrun-by-Naccel matrix where
%      OutputAccel(i,j,k) contains the acceleration recorded for data point
%      i, for run j and accel number k
%
% FreqDomain: A structure containing the frequency response data in the
%    following fields (Here N is the number of identified frequency points
%   - Input Force := An N-by-Nrun-by-Nmeas matrix where InputForce(i,j,k)
%      contains the input force transfered at frequency point i, run j and
%      acceleration number k.
%   - OutputAccel := An N-by-Nrun-by-Nmeas matrix where OutputAccel(i,j,k)
%      contains the output force measured at frequency point i, run j and
%      acceleration number k.
%   - FrequencyHz := An N element array where FrequencyHz(i) contains the
%      frequency at datapoint i.
%   - FreqResponse := An N-by-Nrun-by-Nmeas matrix where FreqResp(i,j,k)
%      contains the frequency response calculated at frequency point i, 
%      run j and acceleration number k.


%% Load experimental information
addpath(dirname)
Info = GVT_Info;
Fs = Info.Fs;
Nmeas = Info.Nmeas;
Nrun = Info.Nrun;
Force_V2lb =  Info.Force_V2lb;
Accel_V2g = Info.Accel_V2g;

%% Load All Data
for i = 1:Nrun
    % Logged Data:
    %   1) timeIn_s is Nt-by-1
    %   2) signalIn is a Nt-by-4 matrix of data
    % Load data and store as:
    %   1) time: Time vector, sec (Nt-by-Nrun)
    %   2) excitation: Excitation input (Nt-by-Nrun)
    %   3) force: Force Transducer output, lbs (Nt-by-Nrun)
    %   4) acccel: Accel output, g's (Nt-by-Nrun-by-Nmeas)
    % Note--This assumes the same time vector for each dataset.
    
    load(fullfile([dirname 'huginn_run' int2str(i) '.mat' ]));
    
    if i==1
        Nt = size(signalIn,1);
        time = timeIn_s;
        excitation = zeros(Nt,Nrun);
        force = zeros(Nt,Nrun);
        accel = zeros(Nt,Nrun,Nmeas);
    elseif size(signalIn,1)~=Nt
        error('# of time points is changing');
    end
    excitation(:,i) = signalIn(:,1);
    force(:,i) = signalIn(:,2)/Force_V2lb;
    accel(:,i,1) = signalIn(:,3)/Accel_V2g;
    accel(:,i,2) = signalIn(:,4)/Accel_V2g;
end

%% Store Time Domain Data
TimeDomain.Time = time;
TimeDomain.InputCommand = excitation;
TimeDomain.InputForce = force;
TimeDomain.OutputAccel = accel;

%% Spectral Analysis
% Perform spectral analysis if user requests this output.
if nargout==3
    for i=1:Nrun
        % Force Input
        u = force(:,i);
        
        % Accel 1
        y = accel(:,i,1);
        [Gtmp1,fHz,Y1,U1] = myspa(y,u,Fs);
        
        % Accel 2
        y = accel(:,i,2);
        [Gtmp2,fHz,Y2,~] = myspa(y,u,Fs);
        
        % Store input/output spectrum and frequency reponse as matrices
        % Note--This assumes the same freq vector for each dataset.
        if i==1
            Nw = numel(fHz);
            U = zeros(Nw,Nrun);
            Y = zeros(Nw,Nrun,Nmeas);
            G = zeros(Nw,Nrun,Nmeas);
        end
        U(:,i)=U1;
        Y(:,i,1)=Y1;
        Y(:,i,2)=Y2;
        G(:,i,1)=Gtmp1;
        G(:,i,2)=Gtmp2;
    end
    FreqDomain.InputForce = U;
    FreqDomain.OutputAccel = Y;
    FreqDomain.FrequencyHz = fHz;
    FreqDomain.FreqResponse = G;
end
end