function u = RK(butcher,a,b,left,rigth,step,s,FV)
syms x t

mesh = [left:step:rigth];
u = zeros(s,length(mesh));
u(:,1) = FV;
omega = zeros(s,s);

for n = 1:length(mesh)-1
    for k = 1:s
        omega(k,:) = RHS(mesh(k)+step*a(k),u(:,k)+step*omega*(butcher(k,:)'));
        omega
    end
  
    u(:,n+1) = u(:,n) + step*omega*(b');
end

end

