function p = efficient_p(U2,U1,U0,r)
%Определяет эффективный порядок по тройке сгущенных сеток
    p = log((U2-U1)/(U1-U0))/log(r);
end

