% Calculate the start points of the trajectory from the terminal state vector xf.
function xs = start_points(xf, M0_max, M0_min, orders, signs, tangents, arctimes)
    n = length(xf);
    N_real = sum(tangents==0);
    xs = zeros([n,N_real+1]);
    xs(:,1) = xf(1:n);
    i_arc = length(tangents);
    for i = 1:N_real
        while (tangents(i_arc)~=0)
            i_arc = i_arc - 1;
        end
        xs(:,i+1) = dynamics_onestep(xs(:,i),cal_u(orders(i_arc),signs(i_arc),M0_max,M0_min),-arctimes(i_arc));
        i_arc = i_arc - 1;
    end
end