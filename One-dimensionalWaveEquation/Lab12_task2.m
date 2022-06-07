%----------Scheme-parametrs---------
global c  tau h a T sigma

c     = 3;
tau   = 0.01;
h     = 6*pi/100;
a     = 6*pi;
T     = 10;
sigma = 1/3;

%----------------Mesh---------------
meshX = [0: h :a];
meshT = [0:tau:T];

%-------------Calculate-------------
U = DifferenceScheme(0);

for i = 1:1:size(meshT,2)
    hold on
    axis([0 meshX(end) -2 2])
    plot(meshX,U(:,i),'r');
    grid on
    pause(1e-5);
    if(i ~= size(meshT,2))
        cla;
    end
    title(['t = ' num2str(meshT(i))]);
    xlabel('x'); ylabel('u(x,t)');  
end