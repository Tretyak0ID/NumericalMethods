function z = mu(x,y)
    z = zeros(length(y),length(x));
    for i=1:length(x)
        for j=1:length(y)
            z(j,i) = -0.01*sin(x(i))+0.05*sin(y(j));
        end
    end
end
