function [ U ] = crosn(k,stPoint,endPoint,step,stTime,endTime,tau,initCond)
gridSpace = stPoint:step:endPoint;
gridTime = stTime:tau:endTime;
U = zeros(length(gridSpace),length(gridTime));
U(1:1:length(gridSpace),1) = initCond(gridSpace);

alpha = (1 + 1i)/2;

matrixDiag = [k/(step^2)*ones(length(gridSpace), 1) -2*k/(step^2)*ones(length(gridSpace), 1) k/step^2*ones(length(gridSpace), 1)];

lymbda = spdiags(matrixDiag, -1:1, length(gridSpace), length(gridSpace));
lymbda(1,:) = 0;
lymbda(1,1) = -2*k/(step^2);
lymbda(1,2) = 2*k/(step^2);
lymbda(end,:) = 0;
lymbda(end,end - 1) = 2*k/(step^2);
lymbda(end,end) = -2*k/(step^2);

E = sparse(eye(length(gridSpace),length(gridSpace)));
matrix = E - alpha * tau * lymbda;
for j = 1:1:(size(gridTime,2))
    w = matrix\(lymbda*U(:,j));
    U(:,j + 1) = U(:,j) + tau * real(w);
end

end