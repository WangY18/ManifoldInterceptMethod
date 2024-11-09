% Solve 3arcs in 3rd-order problems based on Groebner basis.
function [T1,T2,T3] = solution_3arc_3rd_order(x0,xf,orders,signs,M_max,M_min,epsilon)
    u1 = cal_u(orders(1),signs(1),M_max(1),M_min(1));
    u2 = cal_u(orders(2),signs(2),M_max(1),M_min(1));
    u3 = cal_u(orders(3),signs(3),M_max(1),M_min(1)) - u2;
    u2 = u2 - u1;
    x1 = x0(1);
    x2 = x0(2);
    y1 = xf(1) - x0(1);
    y2 = xf(2) - x0(2);
    y3 = xf(3) - x0(3);

    if (abs(u2-u3)>epsilon)
		c6 = u1 * u1 .* (u1 + u2) .* (u1 + u3) .* (u1 + u2 + u3) .* (u1 + u2 + u3);
		c5 = 6 * u1 .* (u1 + u2) .* (u1 + u3) .* (u2 * x1 + u3 * x1 - u1 * y1) .* (u1 + u2 + u3);
		c4 = (15 * y1 * y1 - 6 * u2 * y2 - 6 * u3 * y2) * u1 * u1 * u1 * u1 + (18 * u2 * y1 * y1 + 18 * u3 * y1 * y1 - 12 * u2 * u2 * y2 - 12 * u3 * u3 * y2 + 24 * u2 * u3 * x2 - 12 * u2 * u3 * y2 - 24 * u2 * x1 * y1 - 24 * u3 * x1 * y1) * u1 * u1 * u1 + (12 * u2 * u2 * x1 * x1 + 12 * u3 * u3 * x1 * x1 + 3 * u2 * u2 * y1 * y1 + 3 * u3 * u3 * y1 * y1 - 6 * u2 * u2 * u2 * y2 - 6 * u3 * u3 * u3 * y2 + 12 * u2 * u3 * x1 * x1 + 36 * u2 * u3 * u3 * x2 + 36 * u2 * u2 * u3 * x2 + 15 * u2 * u3 * y1 * y1 - 6 * u2 * u3 * u3 * y2 - 6 * u2 * u2 * u3 * y2 - 24 * u2 * u2 * x1 * y1 - 24 * u3 * u3 * x1 * y1 - 60 * u2 * u3 * x1 * y1) * u1 * u1 + (12 * x2 * u2 * u2 * u2 * u3 + 12 * u2 * u2 * u2 * x1 * x1 + 24 * x2 * u2 * u2 * u3 * u3 + 24 * u2 * u2 * u3 * x1 * x1 - 30 * y1 * u2 * u2 * u3 * x1 + 12 * x2 * u2 * u3 * u3 * u3 + 24 * u2 * u3 * u3 * x1 * x1 - 30 * y1 * u2 * u3 * u3 * x1 + 12 * u3 * u3 * u3 * x1 * x1) * u1 + 9 * u2 * u2 * u2 * u3 * x1 * x1 + 18 * u2 * u2 * u3 * u3 * x1 * x1 + 9 * u2 * u3 * u3 * u3 * x1 * x1;
		c3 = (24 * u2 * y1 * y2 - 24 * u2 * u3 * y3 - 20 * y1 * y1 * y1 + 24 * u3 * y1 * y2) * u1 * u1 * u1 + (36 * u2 * x1 * y1 * y1 - 12 * u3 * y1 * y1 * y1 - 36 * u2 * u3 * u3 * y3 - 36 * u2 * u2 * u3 * y3 - 12 * u2 * y1 * y1 * y1 + 36 * u3 * x1 * y1 * y1 - 24 * u2 * u2 * x1 * y2 - 24 * u3 * u3 * x1 * y2 + 24 * u2 * u2 * y1 * y2 + 24 * u3 * u3 * y1 * y2 - 24 * u2 * u3 * x1 * y2 - 72 * u2 * u3 * x2 * y1 + 24 * u2 * u3 * y1 * y2) * u1 * u1 + (12 * u2 * u2 * x1 * y1 * y1 - 12 * u2 * u3 * u3 * u3 * y3 - 12 * u2 * u2 * u2 * u3 * y3 - 24 * u2 * u2 * u2 * x1 * y2 - 24 * u3 * u3 * u3 * x1 * y2 - 24 * u2 * u2 * u3 * u3 * y3 - 4 * u2 * u3 * y1 * y1 * y1 - 24 * u2 * u2 * x1 * x1 * y1 + 12 * u3 * u3 * x1 * y1 * y1 - 24 * u3 * u3 * x1 * x1 * y1 + 72 * u2 * u3 * u3 * x1 * x2 + 72 * u2 * u2 * u3 * x1 * x2 + 48 * u2 * u3 * x1 * y1 * y1 - 24 * u2 * u3 * x1 * x1 * y1 - 12 * u2 * u3 * u3 * x1 * y2 - 36 * u2 * u3 * u3 * x2 * y1 - 12 * u2 * u2 * u3 * x1 * y2 - 36 * u2 * u2 * u3 * x2 * y1 + 12 * u2 * u3 * u3 * y1 * y2 + 12 * u2 * u2 * u3 * y1 * y2) * u1 + 36 * x2 * u2 * u2 * u2 * u3 * x1 + 8 * u2 * u2 * u2 * x1 * x1 * x1 + 72 * x2 * u2 * u2 * u3 * u3 * x1 - 8 * u2 * u2 * u3 * x1 * x1 * x1 - 36 * y1 * u2 * u2 * u3 * x1 * x1 + 36 * x2 * u2 * u3 * u3 * u3 * x1 - 8 * u2 * u3 * u3 * x1 * x1 * x1 - 36 * y1 * u2 * u3 * u3 * x1 * x1 + 8 * u3 * u3 * u3 * x1 * x1 * x1;
		c2 = (12 * u2 * u2 * y2 * y2 + 72 * y3 * u2 * u3 * y1 + 12 * u2 * u3 * y2 * y2 - 36 * u2 * y1 * y1 * y2 + 12 * u3 * u3 * y2 * y2 - 36 * u3 * y1 * y1 * y2 + 15 * y1 * y1 * y1 * y1) * u1 * u1 + (12 * u2 * u2 * u2 * y2 * y2 + 36 * y3 * u2 * u2 * u3 * y1 - 12 * u2 * u2 * u3 * y2 * y2 - 72 * x2 * u2 * u2 * u3 * y2 - 72 * x1 * y3 * u2 * u2 * u3 - 12 * u2 * u2 * y1 * y1 * y2 + 48 * x1 * u2 * u2 * y1 * y2 + 36 * y3 * u2 * u3 * u3 * y1 - 12 * u2 * u3 * u3 * y2 * y2 - 72 * x2 * u2 * u3 * u3 * y2 - 72 * x1 * y3 * u2 * u3 * u3 - 12 * u2 * u3 * y1 * y1 * y2 + 72 * x2 * u2 * u3 * y1 * y1 + 48 * x1 * u2 * u3 * y1 * y2 + 3 * u2 * y1 * y1 * y1 * y1 - 24 * x1 * u2 * y1 * y1 * y1 + 12 * u3 * u3 * u3 * y2 * y2 - 12 * u3 * u3 * y1 * y1 * y2 + 48 * x1 * u3 * u3 * y1 * y2 + 3 * u3 * y1 * y1 * y1 * y1 - 24 * x1 * u3 * y1 * y1 * y1) * u1 - 36 * y3 * u2 * u2 * u2 * u3 * x1 + 36 * u2 * u2 * u2 * u3 * x2 * x2 - 24 * y2 * u2 * u2 * u2 * x1 * x1 - 72 * y3 * u2 * u2 * u3 * u3 * x1 + 72 * u2 * u2 * u3 * u3 * x2 * x2 + 24 * y2 * u2 * u2 * u3 * x1 * x1 - 72 * u2 * u2 * u3 * x1 * x2 * y1 + 36 * y2 * u2 * u2 * u3 * x1 * y1 + 12 * u2 * u2 * x1 * x1 * y1 * y1 - 36 * y3 * u2 * u3 * u3 * u3 * x1 + 36 * u2 * u3 * u3 * u3 * x2 * x2 + 24 * y2 * u2 * u3 * u3 * x1 * x1 - 72 * u2 * u3 * u3 * x1 * x2 * y1 + 36 * y2 * u2 * u3 * u3 * x1 * y1 + 12 * u2 * u3 * x1 * x1 * y1 * y1 - 12 * u2 * u3 * x1 * y1 * y1 * y1 - 24 * y2 * u3 * u3 * u3 * x1 * x1 + 12 * u3 * u3 * x1 * x1 * y1 * y1;
		c1 = (72 * y3 * u2 * u2 * u3 * y2 - 24 * u2 * u2 * y1 * y2 * y2 + 72 * y3 * u2 * u3 * u3 * y2 - 72 * y3 * u2 * u3 * y1 * y1 - 24 * u2 * u3 * y1 * y2 * y2 + 24 * u2 * y1 * y1 * y1 * y2 - 24 * u3 * u3 * y1 * y2 * y2 + 24 * u3 * y1 * y1 * y1 * y2 - 6 * y1 * y1 * y1 * y1 * y1) * u1 - 72 * x2 * y3 * u2 * u2 * u2 * u3 + 24 * x1 * u2 * u2 * u2 * y2 * y2 - 144 * x2 * y3 * u2 * u2 * u3 * u3 + 72 * x2 * u2 * u2 * u3 * y1 * y2 + 72 * x1 * y3 * u2 * u2 * u3 * y1 - 24 * x1 * u2 * u2 * u3 * y2 * y2 - 24 * x1 * u2 * u2 * y1 * y1 * y2 - 72 * x2 * y3 * u2 * u3 * u3 * u3 + 72 * x2 * u2 * u3 * u3 * y1 * y2 + 72 * x1 * y3 * u2 * u3 * u3 * y1 - 24 * x1 * u2 * u3 * u3 * y2 * y2 - 24 * x2 * u2 * u3 * y1 * y1 * y1 - 24 * x1 * u2 * u3 * y1 * y1 * y2 + 6 * x1 * u2 * y1 * y1 * y1 * y1 + 24 * x1 * u3 * u3 * u3 * y2 * y2 - 24 * x1 * u3 * u3 * y1 * y1 * y2 + 6 * x1 * u3 * y1 * y1 * y1 * y1;
		c0 = 36 * u2 * u2 * u2 * u3 * y3 * y3 - 8 * u2 * u2 * u2 * y2 * y2 * y2 + 72 * u2 * u2 * u3 * u3 * y3 * y3 - 72 * u2 * u2 * u3 * y1 * y2 * y3 + 8 * u2 * u2 * u3 * y2 * y2 * y2 + 12 * u2 * u2 * y1 * y1 * y2 * y2 + 36 * u2 * u3 * u3 * u3 * y3 * y3 - 72 * u2 * u3 * u3 * y1 * y2 * y3 + 8 * u2 * u3 * u3 * y2 * y2 * y2 + 24 * u2 * u3 * y1 * y1 * y1 * y3 + 12 * u2 * u3 * y1 * y1 * y2 * y2 - 6 * u2 * y1 * y1 * y1 * y1 * y2 - 8 * u3 * u3 * u3 * y2 * y2 * y2 + 12 * u3 * u3 * y1 * y1 * y2 * y2 - 6 * u3 * y1 * y1 * y1 * y1 * y2 + y1 * y1 * y1 * y1 * y1 * y1;
        T1 = roots([c6,c5,c4,c3,c2,c1,c0]);
        T1 = T1(abs(imag(T1))<epsilon & real(T1)>0);
        T2_under = -u2 .* (-T1 .* T1 .* T1 * u1 * u1 * u1 * u2 + 3 * T1 .* T1 .* T1 * u1 * u1 * u1 * u3 - T1 .* T1 .* T1 * u1 * u1 * u2 * u2 + 3 * T1 .* T1 .* T1 * u1 * u1 * u2 * u3 + 4 * T1 .* T1 .* T1 * u1 * u1 * u3 * u3 + T1 .* T1 .* T1 * u1 * u2 * u2 * u3 + 2 * T1 .* T1 .* T1 * u1 * u2 * u3 * u3 + T1 .* T1 .* T1 * u1 * u3 * u3 * u3 + 3 * T1 .* T1 * u1 * u1 * u2 * y1 - 9 * T1 .* T1 * u1 * u1 * u3 * y1 + T1 .* T1 * u1 * u2 * u2 * y1 - 2 * x1 * T1 .* T1 * u1 * u2 * u2 - 3 * T1 .* T1 * u1 * u2 * u3 * y1 + 6 * x1 * T1 .* T1 * u1 * u2 * u3 - 4 * T1 .* T1 * u1 * u3 * u3 * y1 + 8 * x1 * T1 .* T1 * u1 * u3 * u3 + 3 * x1 * T1 .* T1 * u2 * u2 * u3 + 6 * x1 * T1 .* T1 * u2 * u3 * u3 + 3 * x1 * T1 .* T1 * u3 * u3 * u3 + 2 * y2 * T1 * u1 * u2 * u2 - 6 * y2 * T1 * u1 * u2 * u3 - 3 * T1 * u1 * u2 * y1 * y1 - 8 * y2 * T1 * u1 * u3 * u3 + 9 * T1 * u1 * u3 * y1 * y1 + 6 * x2 * T1 * u2 * u2 * u3 + 2 * x1 * T1 * u2 * u2 * y1 + 12 * x2 * T1 * u2 * u3 * u3 - 6 * x1 * T1 * u2 * u3 * y1 + 6 * x2 * T1 * u3 * u3 * u3 - 8 * x1 * T1 * u3 * u3 * y1 - 6 * y3 * u2 * u2 * u3 - 2 * y2 * u2 * u2 * y1 - 12 * y3 * u2 * u3 * u3 + 6 * y2 * u2 * u3 * y1 + u2 * y1 * y1 * y1 - 6 * y3 * u3 * u3 * u3 + 8 * y2 * u3 * u3 * y1 - 3 * u3 * y1 * y1 * y1);
        T2_over = (u2 - u3) .* (T1 .* T1 * u1 * u1 + u3 * T1 .* T1 * u1 - 2 * T1 * u1 * y1 + 2 * u3 * x1 * T1 + y1 * y1 - 2 * u3 * y2) .* (y1 * y1 - 2 * u3 * y2 - 2 * u2 * y2 + T1 .* T1 * u1 * u1 + 2 * T1 * u2 * x1 + 2 * T1 * u3 * x1 - 2 * T1 * u1 * y1 + T1 .* T1 * u1 * u2 + T1 .* T1 * u1 * u3);
        % T2_under * T2 + T2_over == 0
        T2 = -T2_over ./ T2_under;
    else
		c3 = u1 * (u1 * u1 + 3 * u1 * u2 + 2 * u2 * u2);
		c2 = (-3 * y1) * u1 * u1 + (6 * u2 * x1 - 3 * u2 * y1) * u1 + 6 * u2 * u2 * x1;
		c1 = (3 * y1 * y1 - 6 * u2 * y2) * u1 + 12 * x2 * u2 * u2 - 6 * x1 * y1 * u2;
		c0 = -12 * y3 * u2 * u2 + 6 * y2 * u2 * y1 - y1 * y1 * y1;
        T1 = roots([c3,c2,c1,c0]);
        T1 = T1(abs(imag(T1))<epsilon & real(T1)>0);
        D2 = 2 * u2 * u2;
        D1 = 2 * T1 * u1 * u2 - 2 * u2 * y1;
        D0 = T1 .* T1 * u1 * u1 + u2 * T1 .* T1 * u1 - 2 * T1 * u1 * y1 + 2 * u2 * x1 * T1 + y1 * y1 - 2 * u2 * y2;
        Delta = D1 .* D1 - 4 * D2 .* D0;
        id = Delta >= 0;
        T1 = T1(id);
        D1 = D1(id);
        D2 = D2(id);
        T2 = (-D1 + sqrt(Delta)) ./ (2 * D2);
    end
    T3 = (u1 * T1 + u2 * T2 - y1) / (-u3);
    id = 0 <= T3 & T3 <= T2 & T2 <= T1;
    T1 = T1(id);
    T2 = T2(id);
    T3 = T3(id);
end