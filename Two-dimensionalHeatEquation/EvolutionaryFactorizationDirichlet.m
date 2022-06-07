function U = EvolutionaryFactorizationDirichlet(uSC_0, uBC_x_0, uBC_x_b, uBC_y_0, uBC_y_a)
global a b T hx hy tau k

meshX = [0: hx:a];
meshY = [0: hy:b];
meshT = [0:tau:T];

Nx = length(meshX);
Ny = length(meshY);
Nt = length(meshT);

U = zeros(Ny, Nx, Nt);

%--------SC-and-BC--------
U(:,:,1)   = uSC_0;
U(1,:,:)   = repmat(uBC_x_0', 1, Nt);
U(end,:,:) = repmat(uBC_x_b', 1, Nt);
U(:,1,:)   = repmat(uBC_y_0,  1, Nt);
U(:,end,:) = repmat(uBC_y_a,  1, Nt);

%-------Differential-operator-------

LambdaX = zeros(Nx);
LambdaY = zeros(Ny);

for i = 2:1:Nx-1
    LambdaX(i,i) = -2*k/hx^2;
    if(i<Nx)
        LambdaX(i+1,i) = 1*k/hx^2;
    end
    if(i>1)
        LambdaX(i-1,i) = 1*k/hx^2;
    end
end

for i = 2:1:Ny-1
    LambdaY(i,i) = -2*k/hy^2;
    if(i<Ny)
        LambdaY(i,i+1) = 1*k/hy^2;
    end
    if(i>1)
        LambdaY(i,i-1) = 1*k/hy^2;
    end
end

Px = (eye(Nx) - (tau/2).*LambdaX);
Py = (eye(Ny) - (tau/2).*LambdaY);

for i = 1:1:Nt-1
    v       = (LambdaY*U(:,:,i)+U(:,:,i)*LambdaX)/Px;
    delta_u = Py\v;
    U(2:end-1, 2:end-1,i+1) = U(2:end-1, 2:end-1,i) + tau*delta_u(2:end-1, 2:end-1);
end
end

