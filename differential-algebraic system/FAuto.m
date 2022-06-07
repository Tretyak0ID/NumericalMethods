function y = FAuto(u)
    global r0 r1 r2 r3 r4 r5 ub;
    y = [
            u(1)/r0 - ue(u(end))/r0; 
            0.01*ff(u(2)-u(3))-ub/r2+u(2)*(1/r1+1/r2);
            u(3)/r3 - ff(u(2)-u(3)); 
            0.99*ff(u(2)-u(3))-ub/r4+u(4)/r4; u(5)/r5;
            1
        ];
end

