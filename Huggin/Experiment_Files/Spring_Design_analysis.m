%% Spring_Design_Analysis.m
% This file contains the code of how the length of the spring was chosen.
% This has no bearing on the GVT experimental analysis.
%
% Note that a separate code exists which contains the measured spring
% stiffness and the rigid body modal frequencies.

%% Clear workspace and conversion constants
clear; close all; clc;

lb2kg = 0.453592;
in2m = 0.0254;
g = 9.81;                                           % in m/s^2
%% Minimum required mass to achieve deflection the full spring
M_min = 1.5;               % in Kg

%% Experiment to obtain stiffness of the full spring
M_stiffness = 2.472;                          % in Kg
T_stiffness = [26.32 26.12 26.20]/20;         % in seconds

T_stiff_avg = mean(T_stiffness);

K_spring = 4*pi^2*M_stiffness/T_stiff_avg^2   % in N/m

%% Analysis for feasibility
% This section contains test of how the spring should be cut
m_aircraft = 40 * lb2kg;
l_spring = 36 * in2m;
h_plate = 38 * in2m;        % Height of the plate from the ground level

% Deflection if full length springs are used:
l_spring_ratio = 1;
d1 = m_aircraft*g/3/(K_spring/l_spring_ratio);
l_req_setup1 = l_spring*l_spring_ratio + d1

% Deflection if 0.5 length springs are used:
l_spring_ratio = 0.5;
d2 = m_aircraft*g/3/(K_spring/l_spring_ratio);
l_req_setup2 = l_spring*l_spring_ratio + d2

% Deflection if 0.3 length springs are used:
l_spring_ratio = 0.3;
d3 = m_aircraft*g/3/(K_spring/l_spring_ratio);
l_req_setup3 = l_spring*l_spring_ratio + d3

% Deflection if 0.25 length springs are used:
l_spring_ratio = 0.25;
d4 = m_aircraft*g/3/(K_spring/l_spring_ratio);
l_req_setup4 = l_spring*l_spring_ratio + d4

%%%%% Note thatl_spring_ration of 0.3 was chosen.