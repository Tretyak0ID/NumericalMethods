function f = F(y,v)
global g d m i v_15 P_0 T_0 w z_t c_d m_3 T_max epsilon T_v S

a = 340.294*sqrt(T(y)/288.15);
f = i*S*(rho(y)*v^2/2)*cx(v/a);

end

