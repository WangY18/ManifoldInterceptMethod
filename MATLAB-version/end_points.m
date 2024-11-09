% Calculate the end points of the trajectory from the initial state vector x0.
function xs = end_points(x0, M0_max, M0_min, orders, signs, tangents, arctimes)
    n = length(x0);
    N_real = sum(tangents==0);
    xs = zeros([n,N_real+1]);
    xs(:,1) = x0;
    i_arc = 1;
    for i = 1:N_real
        while (tangents(i_arc)~=0)
            i_arc = i_arc + 1;
        end
        xs(:,i+1) = dynamics_onestep(xs(:,i),cal_u(orders(i_arc),signs(i_arc),M0_max,M0_min),arctimes(i_arc));
        i_arc = i_arc + 1;
    end
end