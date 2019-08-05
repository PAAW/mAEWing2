%% runSvanbergSQP.m script to run SQP for Svanberg's 5-segment beam. 
% AOE 5064 Homework #3, Problem 5a: SQP for Svanberg beam
Nsegments=5;
Xinitial=5;
X0=repmat(Xinitial,Nsegments,1); vlb=zeros(Nsegments,1);
options=optimset('fmincon');
options.Display='iter';
options.MaxIter=30;
options.TolX=.1;
options.TolFun=.001;
options.TolCon=.001;
[x,output]=sqp(@fSvanbergBeam,X0,options,vlb,[],@gSvanbergBeam);
disp('Final Design Variables, X')
disp(x)