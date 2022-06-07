global g d m i v_15 P_0 T_0 w z_t c_d m_3 T_max epsilon T_v S v_0

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
alpha1 = 0.00001*pi;
alpha2 = pi/4-0.001*pi;

while (1)
    alpha = (alpha1+alpha2)/2;
    [x O] = ode45(@(x,O) ODEsystem(x,O), [0 L], [0 tan(alpha) v_15*cos(alpha) 0]);
    
    if(abs(alpha2-alpha1)<epsilon | O(end,1)==0)
        break
    elseif (O(end,1)<0)
        alpha1 = alpha;
    else
        alpha2 = alpha;
    end
end

hold on; grid on; title('Traectory'); xlabel('x'); ylabel('y');
plot(x, O(:,1))
disp(['angle = ', num2str(alpha)]);

%-------------------------------------------------------------------------
timeMax = O(end - 1,4) - O(end-1,1)*((O(end,4)-O(end-1,4))/(O(end,1)-O(end-1,1)));
vW = 10; 
v_0 = sqrt(v_0^2 + vW^2);
gamma0=tan(alpha);
u0 = v_0*cos(alpha);

[t, O] = ode45(@(t,O) ODEsystem_2(t, O), [0 timeMax], [0,0, gamma0, u0, 0, 0]);
x = O(:,1);
y = O(:,2);
z = O(:,5);

figure
plot3(z,x,y)
grid on
hold on

% отклонение из-за деривации
derDeviat = z(end);

% отклонение из-за ветра
t(end)
windDeviat = vW*t(end)
fprintf('для T = %d \n', T_0);

% отклонение
deviat = windDeviat - derDeviat - 2000*vW/v_0;
fprintf('Отклонение %.2f \n', deviat)

% поправка
corr = (atan(deviat/2000))*3000/pi;
fprintf('горизонтальная угловая поправка %.2f \n', corr)