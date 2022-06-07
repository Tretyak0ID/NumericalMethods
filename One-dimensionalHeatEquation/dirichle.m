figure 

a = 20;
T = 10;
k = 2;

stPoint = 0;
endPoint = a;
step = 0.01;

stTime = 0;
endTime = T;
tau = 0.05;

gridSpace = 0:step:a;
gridTime = 0:tau:T;
timeStep = 0.1;

initCond = @(x)(exp(-((x - 5).^4)) + 0.01*x);

boundCondD1 = @(t)(initCond(0));  
boundCondD2 = @(t)(initCond(a));

[UDirichl] = cros(k,gridSpace(1),gridSpace(end),step,gridTime(1),gridTime(end),tau,boundCondD1,boundCondD2,initCond);
plotAnimation(a,gridSpace,gridTime,UDirichl,timeStep);