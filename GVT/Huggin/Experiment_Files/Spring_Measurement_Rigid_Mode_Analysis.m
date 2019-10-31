%% Spring_Measurement_Rigid_Mode_Analysis:
% This code contains the measurment of the spring stiffness and the 'rigid
% body' modes calculation

%% Run sectionwise only
error('Run the code sectionwise only')

%% Clear Workspace
clear
close all
clc

%% Individual Spring measurement tests
% This section contains the measurement of the stifness of the individual
% springss.

M_Long  = 3.464;       % Mass used for oscillation test

% Spring 1
Time_30cycle_a = [23.08 23.02 23.00];   % Measured time for 30 periods     
Ta = mean(Time_30cycle_a)/30;           % Mean time period
wa = 2*pi/Ta;                           % Frequency (rad/s)
ka = M_Long *wa^2;                     % Spring stiffness


% Spring 2
Time_30cycle_b = [25.82 25.91 25.92];   % Measured time for 30 periods     
Tb = mean(Time_30cycle_b)/30;           % Mean time period
wb = 2*pi/Tb;                           % Frequency (rad/s)
kb = M_Long*wb^2;                      % Spring stiffness

% Spring 3
Time_30cycle_c = [25.42 25.36 25.33 ];  % Measured time for 30 periods     
Tc = mean(Time_30cycle_c)/30;           % Mean time period
wc = 2*pi/Tc;                           % Frequency (rad/s)
kc = M_Long*wc^2   ;                   % Spring stiffness

%% Simplified Rigid Moded Frequency Calculation

in2m = 0.0254;

xa = 21.5 * in2m;                       % Position of Spring1
xb = 36 * in2m;                         % Position of Spring2
xc = 36 * in2m;                         % Position of Spring3
xcg = 30.75*in2m;                       % Position of cg

% Arm length (positive distance)
la = xcg - xa;
lb = xb - xcg;
lc = xc - xcg;

%%%%%%% Longitudinal Model %%%%%%%%%%%%%

mass_AC = 17.45;                        % Mass of the aircraft
Iyy = 1.42;                             % Pitch Inertia of the aircrft

% Longitudinal pitch-plunge motion equaiton
M_Long = [Iyy 0; 0 mass_AC;];           
K_Long = [ka*la + kb*lb + kc*lc   -ka+kb+kc;
     -ka*la + kb*lb + kc*lc   ka+kb+kc;
     ];
 
[v_long,d_long] = eig(K_Long,M_Long);

freq_long = sqrt(d_long)/2/pi


%%%%%%%%%%%%% Lateral Model %%%%%%%%%%%%%%%%%%%%%
% Note that in the lateral case, as the CG lies in the symmetry axis, the
% roll and the plunge motion are decoupled

Ixx = 10.94;                            % Roll Inertia of the aircraft

% Spring Arm Lengths
yb = 8.25 * in2m;
yc = 8.25 * in2m;

% Roll model equations
M_Lat = Ixx;
K_Lat = (kb*yb + kc*yc);

[v_lat,d_lat] = eig(K_Lat,M_Lat);
freq_lat = sqrt(d_lat)/2/pi


%% Measurement of Rigid Body Mode Frequencies 

clear
close all
clc

in2m = 0.0254;
mass_AC = 17.45;                        % Mass of the aircraft

l_ext = 17.5 * in2m; %%
%%%%%%%%%%%%%l_stinger = 8.5 in;
%%%%%% dia = 4.63 mm %%%%%%%%%%

%%%%%%%%%%  Roll %%%%%%%%%%%%%%%%%

I_roll = 10.94;       % Roll Inertia

Time_5cycle_roll = [21.32 21.41 21.37];   % Measured time for 5 periods     
T_roll = mean(Time_5cycle_roll)/5;           % Mean time period
w_roll = 2*pi/T_roll;                         % Frequency (rad/s)
K_roll = I_roll *w_roll^2                     % Stiffness


% Yaw
I_yaw = 13.30;

Time_5cycle_yaw = [25.16 24.84 25.49];   % Measured time for 30 periods     
T_yaw = mean(Time_5cycle_yaw)/5;           % Mean time period
w_yaw = 2*pi/T_yaw;                           % Frequency (rad/s)
K_yaw = I_yaw *w_yaw^2                     % Spring stiffness

% Pitch
I_pitch = 1.42;
Time_5cycle_pitch = [6.38, 6.10 6.18];   % Measured time for 30 periods     
T_pitch = mean(Time_5cycle_pitch)/5;           % Mean time period
w_pitch = 2*pi/T_pitch;                           % Frequency (rad/s)
K_pitch = I_pitch *w_pitch^2                     % Spring stiffness

% Pitch pendulum
Time_5cycle_pitch_pendulum = [7.78 7.81 7.84];   % Measured time for 30 periods        % Measured time for 30 periods     
T_pitch_pendulum = mean(Time_5cycle_pitch_pendulum)/5;           % Mean time period
w_pitch_pendulum = 2*pi/T_pitch_pendulum;                           % Frequency (rad/s)
K_pitch_pendulum = mass_AC *w_pitch_pendulum^2                     % Spring stiffness

% roll pendulum
Time_5cycle_roll_pendulum = [7.15 6.94 6.67];   % Measured time for 30 periods        % Measured time for 30 periods     
T_roll_pendulum = mean(Time_5cycle_roll_pendulum)/5;           % Mean time period
w_roll_pendulum = 2*pi/T_roll_pendulum;                           % Frequency (rad/s)
K_roll_pendulum = mass_AC *w_roll_pendulum^2                     % Spring stiffness

% Heave
Time_5cycle_heave = [5.46 5.5 5.36];   % Measured time for 30 periods        % Measured time for 30 periods     
T_heave = mean(Time_5cycle_heave)/5;           % Mean time period
w_heave = 2*pi/T_heave;                           % Frequency (rad/s)
K_heave = mass_AC *w_heave^2                     % Spring stiffness

