function [ModeShapes] = GVT_ObtainModeShapes(Modal,Info,FreqDomain)
%% Obtains the mode shapes for given modal parameters from the GVT data 
% This function obtains the mode shapes from the frequency domain data for
% a given set of modal frequencies and mode shapes
%
% Inputs:
% Modal: A structure containing modal data
%   - Nmode := Number of mode shapes to be identified
%   - w_mode := A 1-by-Nmode vector containting the identified modal 
%       frequencies in rad/s
%   - zeta := A 1-by-Nmode vector containing the modal damping
%       coefficients for Nmodes
%
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
% FreqDomain: A structure containing the frequency response data in the
%    following fields (Here N is the number of identified frequency points
%    and N_win is the number of identified frequency points in the target
%    range of FreqWin_Hz)
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
%   - G_fit := An N_win-by-Nrun-by-Nmeas matrix where G_fit(i,j,k) contains
%      the frequency response data of the fit model for frequency index i,
%      run j and accelerometer number k.
%   - w_winrad := An N_win-by-1 vector where w_winrad(i) contains the a
%      target frequency at index i.
%
%
% Outputs:
%
% ModeShapes: An Naccel-by-Nmode matrix where ModeShapes(i,j) contains the
%   mode shape data for ith accelerometer and jth mode.


%%


w_mode = Modal.w_mode;
zeta = Modal.zeta;
Nmode = Modal.Nmode;      % Number of modes

Naccel = Info.Naccel;       % Total number of accelerometers
Run2Accel = Info.Run2Accel;

G_fit = FreqDomain.G_fit;
w_winrad = FreqDomain.w_winrad;
%% First approximation of the modes
% This analysis obtains the imaginary part of the frequency response to
% serve as the first approximation of the mode shape

% Interpolating GVT data to get data at modal frequencies
modeshape_gvt_approx=zeros(Naccel,Nmode);
for i = 1:Nmode
    for j = 1:Naccel
        [runj,accelj] = find(Run2Accel==j);
        Gimag = imag(G_fit(:,runj,accelj));
        modeshape_gvt_approx(j,i)=interp1(w_winrad,Gimag,w_mode(i));
    end
end

%% Mode Shape: Approximation 2 (More Accurate)
% This method to obtain the mode shapes is taken from the report by
% Stahle and Forlifer: "Ground Vibration Testing of Complex Structures".
% See p. 86-87 of the report for the corresponding equations.

% The matrix used to extract exact mode shape from the approximate one
Pmat1 = zeros(Nmode);
Pmat2 = zeros(Nmode);
for i=1:Nmode
    for j = 1:Nmode
        % Equivalent to Abhineet's orig. code but different from Fig 4 in
        % Ref.
        if i==j
            Pmat1(i,i) = 1;
        else
            wfrac = w_mode(i)/w_mode(j);
            num = (2*zeta(j))^2*wfrac;
            den = (1-wfrac^2)^2+(2*zeta(j)*wfrac)^2;
            Pmat1(i,j) = num/den;
        end

        % Equation in Fig 4 of Ref
        if i==j
            Pmat2(i,i) = 1;
        else
            wfrac = w_mode(i)/w_mode(j);
            num = (2*zeta(j))^2;
            den = (1-wfrac^2)^2+(2*zeta(j))^2;
            Pmat2(i,j) = num/den;
        end
    end
end

% Edit the following line to select one of the two methods
Pmat = Pmat2;

% Compute the actual mode shape
ModeShapes = zeros(Naccel,Nmode);
for i = 1:Naccel
    ModeShapes(i,:) = (Pmat\modeshape_gvt_approx(i,:)')';
end

end