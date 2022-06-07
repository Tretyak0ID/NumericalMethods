%----------Scheme-parametrs----------
global a b T hx hy tau k
a   = 6*pi;
b   = 4*pi;
T   = 10;
hx  = pi/30;
hy  = pi/30;
tau = 0.1;
k   = 0.2;

%----------Mesh----------
meshX = [0: hx:a];
meshY = [0: hy:b];
meshT = [0:tau:T];

%----------SC-and-BC----------
uSC_0   = mu(meshX, meshY);
uBC_x_0 = mu(meshX, 0);
uBC_x_b = mu(meshX, b);
uBC_y_0 = mu(0, meshY);
uBC_y_a = mu(a, meshY);

U = EvolutionaryFactorizationDirichlet(uSC_0, uBC_x_0, uBC_x_b, uBC_y_0, uBC_y_a);

[X, Y] = meshgrid(meshX, meshY);
for i = 1:1:size(meshT,2)
	mesh(Y, X, U(:,:,i));
	axis([0 meshY(1,end) 0 meshX(1,end) -0.10 0.10]);
    grid on
    pause(1e-5);
    if(i ~= size(meshT,2))
        cla;
    end
end