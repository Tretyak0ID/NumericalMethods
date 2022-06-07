function do = ODEsystem(x,o)
global g d m i v_15 P_0 T_0 w z_t c_d m_3 T_max epsilon T_v S
%y, gamma, u, t
do = zeros(4,1);
do(1) = o(2);
do(2) = -g/o(3)^2;
do(3) = -F(o(1),o(3)*sqrt(1+o(2)^2))/(o(3)*m*sqrt(1+o(2)^2));
do(4) = 1/o(3);
end

