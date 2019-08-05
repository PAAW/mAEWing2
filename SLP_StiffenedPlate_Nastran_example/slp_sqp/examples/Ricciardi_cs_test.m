clear; clc
sqpOptions = optimset('fmincon');
sqpOptions.ComplexStep = 'on';
sqpOptions.DerivativeCheck = 'on';
sqpOptions.Display = 'on';
[x,opts,v,H,status]=sqp(@RicciardiObjCon,[.9;-1],sqpOptions,[-1.5;-1.5],[1.5;1.5]);
