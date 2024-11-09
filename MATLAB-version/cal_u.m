% Calculate the control input u for the arc.
function u = cal_u(order, sign, M0_max, M0_min)
    if (order~=0)
        u = 0;
    elseif (sign>0)
        u = M0_max;
    else
        u = M0_min;
    end
end