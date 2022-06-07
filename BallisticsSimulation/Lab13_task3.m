global g d m i v_15 P_0 T_0 w z_t c_d m_3 T_max epsilon T_v S

g = 9.80665;
d = 12.7*10^(-3);
m = 48.3*10^(-3);
i = 1.0629;
v_15 = 820;
P_0 = 750;
T_0 = 15;
w = 0.5;
z_t = 1.35*10^(-3);
c_d = 0.0423;
m_3 = 0.1744;
T_max = 5;
epsilon = 0.01*pi/3000;

T_v = (T_0+273.15)/(1-3/8*12.7/P_0*w);
S = pi*d^2/4;

L   = 2000;
T_z = T_0;
v_0 = v_15*(1+z_t*(T_z-15));

alpha = 0.01:0.01:1.5;

hold on; grid on; title('len'); xlabel('alpha'); ylabel('len');

for j=1:length(alpha)
    [x O] = ode45(@(x,O) ODEsystem(x,O), [0 7000], [0 tan(alpha(j)) v_15*cos(alpha(j)) 0]);
    for k = 1:length(x)-1
        if(O(k,1)*O(k+1,1)<0)
            plot(alpha(j),x(k),'.b');
        end
    end
end
