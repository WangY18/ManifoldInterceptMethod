% Solve tangent markers in 3rd-order problems based on Groebner basis.
% If positive_time is true, then x0->arc1->arc2->(3,2); else, (3,2)->arc2->arc1->x0.
function [T1,T2] = solution_2arc_tangent_3rd_order(x0,orders,signs,M_max,M_min,positive_time,epsilon)
    u1 = cal_u(orders(1),signs(1),M_max(1),M_min(1));
    u2 = cal_u(orders(2),signs(2),M_max(1),M_min(1)) - u1;
    x1 = x0(1);
    x2 = x0(2);
    y2 = - x0(2);
    if (positive_time==signs(2))
        y3 = M_max(4) - x0(3);
    else
        y3 = M_min(4) - x0(3);
    end
        
	c6 = u1 * u1 * (u1 + u2);
	c5 = 6 * u1 * x1 * (u1 + u2);
	c4 = (-6 * y2) * u1 * u1 + (12 * x1 * x1 + 12 * u2 * x2) * u1 + 9 * u2 * x1 * x1;
	c3 = (-12 * u2 * y3 - 24 * x1 * y2) * u1 + 8 * x1 * x1 * x1 + 36 * u2 * x2 * x1;
	c2 = -24 * x1 * x1 * y2 - 36 * u2 * y3 * x1 + 36 * u2 * x2 * x2 + 12 * u1 * y2 * y2;
	c1 = 24 * x1 * y2 * y2 - 72 * u2 * x2 * y3;
	c0 = -8 * y2 * y2 * y2 + 36 * u2 * y3 * y3;
    T1 = roots([c6,c5,c4,c3,c2,c1,c0]);
    T1 = T1(abs(imag(T1))<epsilon & real(T1)>0);
    T2 = (2 * y2 - u1 * T1 .* T1 - 2 * x1 * T1) / u2;
    id = T2 >= 0;
    T1 = T1(id);
    T2 = T2(id);
    if positive_time
        T2 = sqrt(T2);
    else
        T2 = -sqrt(T2);
    end
    id = T2 <= T1;
    T1 = T1(id);
    T2 = T2(id);
    xf1 = x1 + u1 * T1 + u2 * T2;
    if (positive_time)
        id = xf1 >= 0;
    else
        id = xf1 <= 0;
    end
    T1 = T1(id);
    T2 = T2(id);
end