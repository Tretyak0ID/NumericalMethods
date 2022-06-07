function [mesh,u] = DiffAlgROSAuto(left, rigth, tau, G, alpha, u0)
%Решение дифференциально-алгебраических систем на основании одностадийной
%схемы Розенброка

s = length(u0);
h = 10^(-5);

mesh = [left:tau:rigth];
u    = zeros(s,length(mesh));
u(:,1) = u0;
dFdu = zeros(s);

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
        
        rd                  = FAuto(first_point)-FAuto(second_point);
        dFdu(:,k)           = rd/(2*h*u(k,n));
    end
    
    omega    = (G-alpha*tau*dFdu)\FAuto(u(:,n));
    u(:,n+1) = u(:,n) + tau*real(omega);
    
end

end

