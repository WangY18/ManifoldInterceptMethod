% Plan a 2nd-order trajectory from x0 to xf under constraint.
% Only: 00, 010
% If direction is 1/-1, then directly consider the case where x0 is lower/higher than xf.
function [orders,signs,tangents,arctimes] = plan_2nd_order(x0,xf,M_max,M_min,flag_consider_position,direction,epsilon)
    if (direction==0)
        [orders_try,signs_try,tangents_try,arctimes_try] = plan_1st_order(x0(1),xf(1),M_max(1),M_min(1));
        xs_try = end_points(x0,M_max(1),M_min(1),orders_try,signs_try,tangents_try,arctimes_try);
        if (max(abs(xs_try(:,end)-xf))<epsilon)
            if (flag_consider_position && ~feasible(x0,M_max,M_min,orders_try,signs_try,tangents_try,arctimes_try,epsilon))
                orders = [];
                signs = [];
                tangents = [];
                arctimes = [];
            else
                orders = orders_try;
                signs = signs_try;
                tangents = tangents_try;
                arctimes = arctimes_try;
            end
            return
        end
        if (xs_try(2,end)<xf(2))
            direction = 1;
        else
            direction = -1;
        end
    end

    if (direction<0) % x0 is higher than xf. It must be feasible.
        [orders,signs,tangents,arctimes] = plan_2nd_order(-x0,-xf,-M_min,-M_max,flag_consider_position,1,epsilon);
        signs = ~signs;
        return
    end

    % x0 is lower than xf
    DeltaX2 = xf(2) - x0(2) ...
        + (M_max(2) * M_max(2) - xf(1) * xf(1)) / (2 * M_min(1))...
        - (M_max(2) * M_max(2) - x0(1) * x0(1)) / (2 * M_max(1));
    % Consider 010
    if (~isinf(M_max(2)) && DeltaX2>0)
        orders = [0,1,0];
        signs = [true,true,false];
        tangents = [0,0,0];
        arctimes = [(M_max(2)-x0(1))/M_max(1), DeltaX2/M_max(2), (xf(1)-M_max(2))/M_min(1)];
        if (flag_consider_position && ~feasible(x0,M_max,M_min,orders_try,signs_try,tangents_try,arctimes_try,epsilon))
            orders = [];
            signs = [];
            tangents = [];
            arctimes = [];
        end
        return
    end
    % Consider 00
    T2 = max(roots([M_min(1)*(M_min(1)-M_max(1)),2*xf(1)*(M_max(1)-M_min(1)),xf(1)*xf(1)-x0(1)*x0(1)+2*M_max(1)*(x0(2)-xf(2))]));
    T1 = (xf(1)-x0(1)-(M_min(1)-M_max(1))*T2)/M_max(1); % T1 > T2 > 0 holds.
    orders = [0,0];
    signs = [true,false];
    tangents = [0,0];
    arctimes = [T1-T2,T2];
end