function [U] = DifferenceScheme(flag)
global c  tau h a T sigma

%----------------Mesh-------------------
meshX = [0: h :a];
meshT = [0:tau:T];

Nx    = length(meshX);
Nt    = length(meshT);

U     = zeros(Nx,Nt);

%----------------SC-and-BC---------------
if(flag == 1)
    U(:,1) = mu3_1(meshX);
else
    U(:,1) = mu3_2(meshX);
end
U(1,:)   = mu1(meshT);
U(end,:) = mu2(meshT);

%----------Differential-operator---------
diags = [(c/h)^2*ones(Nx, 1) -2*(c/h)^2*ones(Nx, 1) (c/h)^2*ones(Nx, 1)];
Lambda = spdiags(diags, -1:1, Nx, Nx);
Lambda = tau^2*Lambda;
Lambda(:,1) = 0;
Lambda(:,end) = 0;
E = sparse(eye(Nx));
P = (E - sigma*Lambda)\Lambda;

%-----------First-layer-approx-----------

if(flag == 1)
    u_layer_1 = mu3_1(meshX) + tau*mu4(meshX) + 1/2*((Lambda*(mu3_1(meshX))')' + tau^2*f(meshX,0));
    U(2:end-1,2) =  u_layer_1(2:end-1);
else
    u_layer_2 = mu3_2(meshX) + tau*mu4(meshX) + 1/2*((Lambda*(mu3_2(meshX))')' + tau^2*f(meshX,0));
    U(2:end-1,2) =  u_layer_2(2:end-1);
end

%----------------Caculate----------------
for i = 2:1:Nt-1
    w = P*U(:,i);
    U(2:end-1,i+1) = w(2:end-1) + 2.*U(2:end-1,i) - U(2:end-1,i-1);
end

end

