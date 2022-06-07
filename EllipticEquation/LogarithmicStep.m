function [tau, meshT] = LogarithmicStep(N)
global a b hx hy k epsilon

lambda_x1 = (pi/a)^2*k;
lambda_y1 = (pi/b)^2*k;
lambda_xN = 4*k/hx^2;
lambda_yN = 4*k/hy^2;

lambda_min = min(abs(lambda_x1), abs(lambda_y1));
lambda_max = max(abs(lambda_xN), abs(lambda_yN));

tau_min = 2/lambda_max;
tau_max = 2/lambda_min;

S = round(1/4*log(1/epsilon)*log(2*N/pi));

tau = zeros(1,S);
for s=1:S
    tau(s) = tau_min*(tau_max/tau_min)^((s-1)/(S-1));
end

meshT = zeros(1,S+1);
for s=2:S
    meshT(s) = meshT(s-1)+tau(s);
end
    
end

