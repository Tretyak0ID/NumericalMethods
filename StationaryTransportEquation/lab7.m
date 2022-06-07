clc
clear
close all

x = 0:0.1:100;
t = 0:0.01:1;
c = 50;
u = 1./(1 + ((x - c*t(end-1) - 20)/10).^10);
u1 = neiav(t, x, c);
u2 = complexf(t, x, c);
u3 = sum05(t, x, c);

for i = 1 : length(t)-1
    clf
    hold on
    grid on
    plot(x, u2(i, :), x, u1(i, :),x, u3(i, :), x, u)
    legend('complex','implicit', 'half-sum', 'solve') 
    pause(1e-6)
end

max(abs(u2(end-1, :) - u))