function y =  F5(u)

global N h k0 k1 x ff

y = zeros(N+1,1);
y(1) = u(1);
y(N+1) = u(N+1);

for n = 2:N
    y(n) = (u(n+1) - u(n))*(k0 + k1*(u(n)^2+u(n+1)^2)/2) -(u(n) - u(n-1))*(k0 + k1*(u(n)^2+u(n-1)^2)/2) -h^2*ff(x(n));
end
end

