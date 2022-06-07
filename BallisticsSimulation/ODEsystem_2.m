function do = ODEsystem_2(t, o)
global g d m i v_15 P_0 T_0 w z_t c_d m_3 T_max epsilon T_v S v_0

do = zeros(6,1);
do(1) = o(4);
do(2) = o(3)*o(4);
do(3) = -g/o(4);
do(4) = -F(o(2),o(4)*sqrt(1+o(3)^2))/(m*sqrt(1+o(3)^2));
do(5) = o(6)*o(4)*pi*v_0*c_d;
do(6) = exp(-m_3.*t)/(o(4)^2*(1+o(3)^2));
end

