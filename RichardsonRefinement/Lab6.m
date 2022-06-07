%Параметры задачи
global a b k0 k1 n num h
k0 = 1;
k1 = 0.5;
n = 8;
a = 0;
b = 1;
num = 7;
h = (b - a)/n;
eps = 10^-12;

%Собственно рассчет (рншение + собственные значения)
u0 = zeros(1,n+2);
[x,u,lambda] = Newton(u0, eps);

figure(1)
hold on; grid on
title('Solution')
plot(x,u(1:end-1),'m')

figure(2)
hold on; grid on
title('\lambda')
plot(1:length(lambda), lambda,'.-b')

%Эффективный порядок метода
p = zeros(length(lambda)-2, 1);
for i = 1:length(lambda)-2
    p(i) = -log2(abs((lambda(i+2) - lambda(i+1))/(lambda(i+1) - lambda(i))));
end

figure(3)
hold on; grid on
title('Degree')
plot(p,'.-c')

% уточнение по Ричадрсону
RichardsonRefinement(lambda)