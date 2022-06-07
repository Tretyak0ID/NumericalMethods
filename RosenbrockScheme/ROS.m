function [mesh,u] = ROS(left, rigth, tau, alpha, u0)
s = length(u0);
h = 10^(-5);

mesh = [left:tau:rigth];
u    = zeros(s,length(mesh));
u(:,1) = u0;
dRHSdu = zeros(s);

for n = 1:length(mesh)-1
    
    for k = 1:s
        first_point        = zeros(1,s);
        second_point       = zeros(1,s);
        
        first_point(1:k-1) = u(1:k-1,n);
        first_point(k)     = u(k,n)*(1+h);
        first_point(k+1:s) = u(k+1:s,n);
        
        second_point(1:k-1) = u(1:k-1,n);
        second_point(k)     = u(k,n)*(1-h);
        second_point(k+1:s) = u(k+1:s,n);
        
        rd                  = RHS(mesh(n),first_point)-RHS(mesh(n), second_point);
        dRHSdu(:,k)         = rd/(2*h*u(k,n));
    end

    
    omega    = (eye(s)-alpha*tau*dRHSdu)\RHS(mesh(n)+tau/2,u(:,n));
    u(:,n+1) = u(:,n) + tau*real(omega);
    
end

end

