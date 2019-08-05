function [speed,density,temperature,dynamic_pressure,unit] = atmosphere_property_altitude(alti,Mach)
% this subroutine is used to generate atomosphere properties for NASTRAN,
% use lbs for all weights, and then in NASTRAN, setup WTMASS to convert
% them back to mass as necessary.
%
%
% https://www.fighter-planes.com/jetmach1.htm
% alt - ft
% speed of sound (sos) - knots
%
alt_sos = [0 	 	661 	
 5000 	 	650 	
10000 	 	638 	
15000 	 	626 	
20000 	 	614 	
25000 	 	602 	
30000 	 	589 	
35000 	 	574 	
40000 	 	573 	
45000 	 	573 	
50000 	 	573 	
55000 	 	573 	
60000 	 	573 	];


specified_sos_alt = interp1(alt_sos(:,1),alt_sos(:,2),alti);
knots_to_fps = 1.68781;

speed = specified_sos_alt*Mach*knots_to_fps;



%%  https://www.engineeringtoolbox.com/standard-atmosphere-d_604.html
% Geo-potential Altitude above Sea Level
% - h -
% (ft)	Temperature
% - t -
% (oF)	Acceleration of Gravity
% - g -
% (ft/s2)	Absolute Pressure
% - p -
% (lb/in2)	Density
% - ? -
% (10-4 slugs/ft3)
% Dynamic Viscosity
% - ? -
% (10-7 lb s/ft2)
% 
% ALTITUDE   
standard_US = ...
[0 	     59 	32.174 	14.696 	23.77 	        3.737
5000 	41.17 	32.159 	12.228 	20.48 	        3.637
10000 	23.36 	32.143 	10.108 	17.56 	        3.534
15000 	5.55 	32.128 	8.297 	14.96 	        3.430
20000 	-12.26 	32.112 	6.759 	12.67 	        3.324
25000 	-30.05 	32.097 	5.461 	10.66 	        3.217
30000 	-47.83 	32.082 	4.373 	8.91 	        3.107
35000 	-65.61 	32.066 	3.468 	7.38 	        2.995
40000 	-69.70 	32.051 	2.730 	5.87 	        2.969
45000 	-69.70 	32.036 	2.149 	4.62 	        2.969
50000 	-69.70 	32.020 	1.692 	3.64 	        2.969
60000 	-69.70 	31.990 	1.049 	2.26 	        2.969
70000 	-67.42 	31.959 	0.651 	1.39 	        2.984
80000 	-61.98 	31.929 	0.406 	0.86 	        3.018
90000 	-56.54 	31.897 	0.255 	0.56 	        3.052
100000 	-51.10 	31.868 	0.162 	0.33 	        3.087
150000 	19.40 	31.717 	0.020 	0.037 	        3.511
200000 	-19.78 	31.566 	0.003 	0.0053 	        3.279
250000 	-88.77 	31.415 	0.000 	0.00065 	2.846
];

slugs_to_lbs = 32.174; % 
standard_US(:,5) = standard_US(:,5)*1E-4; standard_US(:,5) = standard_US(:,5)*slugs_to_lbs ; % convert to lbs/ft^3



standard_US(:,6) = standard_US(:,6)*1E-7;


density = interp1(standard_US(:,1),standard_US(:,5),alti);

dynamic_pressure = 0.5*density*speed^2;


temperature = interp1(standard_US(:,1),standard_US(:,2),alti);


unit={'ft/s','lbs/ft^3','psf','F'};
