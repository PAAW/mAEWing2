%% Regression_test script to run slp_trust and sqp examples

%% Log file
clear; clc; close all
delete regression_test.txt
diary  regression_test.txt

%% Run scripts
run2barFox
regression_check( xf, [1.87835663, 20.23690725] );

Ricciardi_cs_test
runBarnes
runBeamGVslp
runHaftka4p2p1
runHaftka6p3p1slp
runRosenSuzuki
runSpringExampleTrust

runSvanbergSLP
xSvanberg = [
   6.01858007
   5.30292693
   4.49935286
   3.50102575
   2.15179983];
regression_check( x, xSvanberg, 0.03 );

runSvanbergSQP
xSvanberg = [
   6.01858007
   5.30292693
   4.49935286
   3.50102575
   2.15179983];
regression_check( x, xSvanberg, 0.06 );

tsuite

diary off