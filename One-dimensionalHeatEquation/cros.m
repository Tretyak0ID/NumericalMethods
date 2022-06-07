function [ U ] = cros(k,stPoint,endPoint,step,stTime,endTime,tau,boundCond1,boundCond2,initCond)

gridSpace = stPoint:step:endPoint;
gridTime = stTime:tau:endTime;

U = zeros(length(gridSpace),length(gridTime));
U(1,1:1:length(gridTime)) = boundCond1(gridTime);
U(end,1:1:length(gridTime)) = boundCond2(gridTime);
U(1:1:length(gridSpace),1) = initCond(gridSpace);

alpha = (1 + 1i)/2;

matrixDiag = [k/(step^2)*ones(length(gridSpace), 1) -2*k/(step^2)*ones(length(gridSpace), 1) k/step^2*ones(length(gridSpace), 1)];
lymbda = spdiags(matrixDiag, -1:1, length(gridSpace), length(gridSpace));
lymbda(1,:) = 0;
lymbda(end,:) = 0;

E = eye(length(gridSpace),length(gridSpace));

for j = 1:1:(size(gridTime,2))
    matrix = E - alpha * tau * lymbda;
    w = matrix\(lymbda*U(:,j));
    U(2:1:(end-1),j + 1) = U(2:(end-1),j) + tau * real(w(2:(end-1)));
end

end