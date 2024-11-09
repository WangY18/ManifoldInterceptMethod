% Plan an $n$th-order trajectory from x0 to xf under constraint.
% If direction is 1/-1, then directly consider the case where x0 is lower/higher than xf.
function [orders,signs,tangents,arctimes] = plan_nth_order(x0,xf,M_max,M_min,flag_consider_position,direction,epsilon)
    if nargin<7
        epsilon = 0;
    end
    n = length(x0);
    if (n==1)
        [orders,signs,tangents,arctimes] = plan_1st_order(x0,xf,M_max(1),M_min(1));
        return
    elseif (n==2)
        [orders,signs,tangents,arctimes] = plan_2nd_order(x0,xf,M_max,M_min,flag_consider_position,direction,epsilon);
        return
    elseif (n==3)
        [orders,signs,tangents,arctimes] = plan_3rd_order(x0,xf,M_max,M_min,flag_consider_position,direction,epsilon);
        return
    end
    % TODO: higher-order version
    orders = [];
    signs = [];
    tangents = [];
    arctimes = [];
end