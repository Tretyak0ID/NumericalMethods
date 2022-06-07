function u_new_1 = ConservativeScheme(u_old_1, u_new_0)
global h tau
    u_new_1 =sqrt(h^2/tau^2 + 2*h/tau*u_old_1+u_new_0^2)-h/tau;
end

