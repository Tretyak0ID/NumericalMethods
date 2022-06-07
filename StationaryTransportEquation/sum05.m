function [ u ] = sum05( t, x, c)

u = zeros(size(t, 2), size(x, 2));
u(:, 1) = 1./(1 + (-5 * t - 2).^10);
u(1, :) = 1./(1 + ((x - 20)/10).^10);

h = x(2) - x(1);
tau = (t(2) - t(1));

for i = 2 : size(u, 1)
    for j = 2 : size(u, 2)
        u(i, j) = (h*(u(i-1, j) + u(i-1, j-1) - u(i, j-1)) + c*tau*(u(i, j-1) + u(i-1, j-1) - u(i-1, j)))/(h + c * tau);
    end
end

end