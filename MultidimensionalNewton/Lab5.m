global N h k0 k1 x ff
ff = @(x) 100*exp(-(10*(x-0.5)).^2);

%Рассчетная сетка
a = -1;
b = 1;
N = 128;
h = (b-a)/N;
x = a:h:b;

%Tочность
epsilon = 10^(-12);

%Параметры схемы
k0 = 1;
k1 = 0.05;

u = zeros(size(x, 2), 1);

u = Newton(u, epsilon);

hold on; grid on;
plot(x,u, 'm')
title('non-liner scheme solution')