function [U,residual] = EvolutionaryFactorizationDirichlet(uSC_0, uBC_x_0, uBC_x_b, uBC_y_0, uBC_y_a)
global a b hx hy k

%----------------Mesh-------------------
meshX        = [0: hx:a];
meshY        = [0: hy:b];
[tau, meshT] = LogarithmicStep(length(meshX));

Nx = length(meshX);
Ny = length(meshY);
Nt = length(meshT);

U        = zeros(Ny, Nx, Nt);
residual = zeros(1,Nt - 1);

%----------------SC-and-BC---------------
U(:,:,1)   = uSC_0;
U(1,:,:)   = repmat(uBC_x_0', 1, Nt);
U(end,:,:) = repmat(uBC_x_b', 1, Nt);
U(:,1,:)   = repmat(uBC_y_0,  1, Nt);
U(:,end,:) = repmat(uBC_y_a,  1, Nt);

%----------Differential-operator---------
Bx = [k/hx^2*ones(Nx, 1) -2*k/hx^2*ones(Nx, 1) k/hx^2*ones(Nx, 1)];
LambdaX = spdiags(Bx, -1:1, Nx, Nx);
LambdaX(:,1) = 0;
LambdaX(:,end) = 0;

By = [k/hy^2*ones(Ny, 1) -2*k/hy^2*ones(Ny, 1) k/hy^2*ones(Ny, 1)];
LambdaY = spdiags(By, -1:1, Ny, Ny);
LambdaY(1,:) = 0;
LambdaY(end,:) = 0;

Ex = sparse(eye(Nx));
Ey = sparse(eye(Ny));

%----------------Caculate----------------
for i = 1:1:Nt-1
    nev                     = [zeros(1,Nx); U(2:end-1,:,i)*LambdaX; zeros(1,Nx)] + [zeros(Ny,1), LambdaY*U(:,2:end-1,i), zeros(Ny,1)];
    v                       = nev/(Ex - tau(i)/2.*LambdaX);
    delta_u                 = (Ey - tau(i)/2.*LambdaY)\v;
    U(2:end-1, 2:end-1,i+1) = U(2:end-1, 2:end-1,i) + tau(i)*delta_u(2:end-1, 2:end-1);
    residual(i)             = max(max(abs(nev)));
end
end

