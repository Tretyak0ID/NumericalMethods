%----------Scheme-parametrs---------
global a b hx hy k epsilon
a   = 6*pi;
b   = 4*pi;
hx  = pi/30;
hy  = pi/30;
k   = 0.2;
epsilon = 10^(-12);

%----------------Mesh---------------
meshX = [0: hx:a];
meshY = [0: hy:b];

%---------Logarithmic-step----------
[tau, meshT] = LogarithmicStep(length(meshX));

%-------------SC-and-BC-------------
uSC_0   = mu(meshX, meshY);
uBC_x_0 = mu(meshX, 0);
uBC_x_b = mu(meshX, b);
uBC_y_0 = mu(0, meshY);
uBC_y_a = mu(a, meshY);

%-------------Calculate-------------
[U, residual] = EvolutionaryFactorizationDirichlet(uSC_0, uBC_x_0, uBC_x_b, uBC_y_0, uBC_y_a);

hold on; grid on;
[X, Y] = meshgrid(meshX, meshY);
mesh(Y, X, U(:,:,end));

figure; hold on; grid on;
plot(1:length(meshT)-1,log10(residual))
xlabel('Number of iteration')
ylabel('Discrepancy ')