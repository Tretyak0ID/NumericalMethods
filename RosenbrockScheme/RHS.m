function  y = RHS(t,u)
    y = [-50*(u(1)-cos(t))+10*u(2);1.2*u(1)-u(2)*u(1)];
end

