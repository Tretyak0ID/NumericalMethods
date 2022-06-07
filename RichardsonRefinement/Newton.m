function [x,u,lambda] = Newton(u0, eps)
global a b k0 k1 n num h

%Начальные значения для решения
lambda0 = 40;
du = ones(1,n+2);
u0(end) = lambda0;
u0(end/2+1) = h;
u0(end/2-1) = -h;
x = a:h:b;
u = u0;
lambda = zeros(num + 1, 1);
lambda(1) = lambda0; 

for j = 1 : num
    while(max(abs(du)) > eps)
        
        %Формирование новой матрицы Якоби
        dF = zeros(n+2);
               
        dF = dF + diag([(k0 + k1*(u(2:end-2).^2 + u(1:end-3).^2)/2) - (u(2:end-2) - u(1:end-3)).*(k1*u(1:end-3)) 0 0], -1) + diag([1 (-(k0 + k1*(u(2:end-2).^2 + u(3:end-1).^2)/2) + (u(3:end-1) - u(2:end-2)).*(k1*u(2:end-2))-(k0 + k1*(u(2:end-2).^2 + u(1:end-3).^2)/2) - (u(2:end-2) - u(1:end-3)).*(k1*u(2:end-2)) + h^2*u(end)) 1 0], 0) + diag([0 (k0 + k1*(u(2:end-2).^2 + u(3:end-1).^2)/2) + (u(3:end-1) - u(2:end-2)).*(k1*u(3:end-1)) 0], 1); 
        dF(2:end-2, end) = h^2*u(2:end-2);
        dF(end, end/2-1) = -1;
        dF(end, end/2+1) = 1;

        %Метод Ньютона
        F = zeros(n+2, 1);
        F(1) = u(1);    %ГУ
        F(2:end-2) = (u(3:end-1) - u(2:end-2)).*(k0 + k1*(u(3:end-1).^2 + u(2:end-2).^2)/2)-(u(2:end-2) - u(1:end-3)).*(k0 + k1*(u(1:end-3).^2 + u(2:end-2).^2)/2)+ h^2*u(2:end-2)*u(end);
        F(end-1) = u(end-1); %ГУ
        F(end) = u(end/2+1)-u(end/2-1)-2*h; %ГУ на лямбду
        F = -F;

        G = dF\F;
        du = G(:, end);
        u = u + du';
    end
    
    lambda(j+1) = u(end);
    
    %Пересчет значений для вычисления слудующего приближения
    u0 = u;
    n = 2*n;
    h = (b-a)/n;
    x = a:h:b;
    u = zeros(1,n+2);
    du = ones(1,n+2);
    u(1:2:end-1) = u0(1:end-1);
    u(end) = u0(end);
    u(2:2:end-2) = (u0(1:end-2) + u0(2:end-1))/2;
end
end

