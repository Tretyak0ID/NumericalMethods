%----------Scheme-parametrs---------
global c  tau h a T sigma

c     = 3;
a     = 6*pi;
T     = 10;
sigma = 1/3;

%----------------Mesh---------------
meshX = [0: h :a];
meshT = [0:tau:T];


layers = 8;
err = zeros(1,layers);
for i = 0:1:(layers - 1)
    h = (6*pi/16)/(2^i);
    tau = (1/16)/(2^i);
    
    U = DifferenceScheme(0);
    if(i ~= 0)
        err(i) = max(abs(U(1:2:end,end) - U_old(:,end)));
    end
    U_old = U;
end

p = zeros(1,layers - 2);
for i = 1:1:(layers - 2)
    p(i) = -log2((err(i + 1))/(err(i)));
end
figure()
plot(1:(layers - 2),p,'r')
grid on