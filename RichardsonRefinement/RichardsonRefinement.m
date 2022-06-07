function lambda_out = RichardsonRefinement(lambda)
global num

lambda_pr = lambda(2:end); %промежуточное лямбда

for j = 1:floor(num-2)
    p = zeros(length(lambda_pr)-2, 1);
    
    for i = 1:length(lambda_pr)-2
        p(i) = -log2(abs((lambda_pr(i+2) - lambda_pr(i+1))/(lambda_pr(i+1) - lambda_pr(i))));
    end
    
    R = zeros(length(lambda_pr)-1, 1);
    
    for i = 1:length(R)
        R(i) = (lambda_pr(i+1) - lambda_pr(i)) / (2^round(p(end)) - 1);
    end
    lambda_pr = lambda_pr(2:end) + R;
end
lambda_out = (2^round(p(end))*lambda_pr(2) - lambda_pr(1)) / (2^round(p(end)) - 1);
end

