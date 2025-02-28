% Plan a 1st-order trajectory from x0 to xf under constraint. (trivial)
function [orders,signs,tangents,arctimes] = plan_1st_order(x0,xf,M0_max,M0_min)
    orders = 0;
    tangents = 0;
    if (x0 < xf)
        signs = true;
        arctimes = (xf - x0) / M0_max;
    else
        signs = false;
        arctimes = (xf - x0) / M0_min;
    end
end