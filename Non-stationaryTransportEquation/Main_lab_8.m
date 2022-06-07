global x_min x_max h t_min t_max tau

x_min = 0;
x_max = 100;
h     = 1;
t_min = 0;
t_max = 100;
tau   = 0.1;

x = [x_min: h:   x_max];
t = [t_min: tau: t_max];
mesh = zeros(length(x),length(t));

for i=1:1:length(x)
    mesh(i,1) = BC(x(i));
end

for i=1:1:length(t)
    mesh(1,i) = mesh(1,1);
end

for i=2:1:length(t)
    for j=2:1:length(x)
        mesh(j,i) = ConservativeScheme(mesh(j,i-1), mesh(j-1,i));
    end
end

for i=1:length(t)
    clf;
    hold on; grid on;
    plot(x,mesh(:,i))
    pause(1e-6);
end