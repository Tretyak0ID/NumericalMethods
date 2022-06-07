function u = Newton(u, epsilon)
global N
    dx = inf;
    J = zeros(N+1);
    rho = 10^-5;

    while(norm(dx) > epsilon)
        for j = 2:N
            dU = zeros(N+1,1);
            dU(j) = rho;
            J(:,j) = (F5(u + dU) - F5(u - dU))/(2*rho);
        end

        J(1,1) = 1;
        J(N+1,N+1) = 1;
        dx = J\-F5(u);
        u = u + dx;
    end

end

