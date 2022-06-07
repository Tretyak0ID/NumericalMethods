clc; clear; figure
a = 20; %заданные условия
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

initCond = @(x)(exp(-((x - 5).^4)) + 0.01*x); %Начальное условие

[ UNeuman ] = crosn(k,gridSpace(1),gridSpace(end),step,gridTime(1),gridTime(end),tau,initCond);
plotAnimation(a,gridSpace,gridTime,UNeuman,timeStep);