function [] =  plotAnimation(a,gridSpace,gridTime,matrixUSC,timeStep)
for i = 1:1:length(gridTime)
    hold on
    axis([0 a 0 1])
    plot(gridSpace,matrixUSC(:,i),'b');
    grid on
    pause(timeStep);
    if(i ~= length(gridTime))
        cla;
    end
end
end