% Check if the trajectory is feasible.
function flag_feasible = feasible(x0, M_max, M_min, orders, signs, tangents, arctimes, epsilon)
    x0 = x0(:);
    n = length(x0);
    if (any(x0<=M_min(2:n+1)-epsilon | x0>=M_max(2:n+1)+epsilon))
        flag_feasible = false;
        return
    end
    T = ones([n+1,1]); % T(j+1) = dt^j/(j!)
    for i = 1:length(orders)
        if (tangents(i)~=0)
            continue
        end
        % Check the local minimum or maximum of x(k+1) along the trajectory.
        u = cal_u(orders(i),signs(i),M_max(1),M_min(1));
        for k = 1:n-1
            T(k+1) = T(k)/k;
            % Find the root of x(k): $k$th-order polynomial
            rs = roots([T(k+1)*u;T(k:-1:1).*x0(1:k)]);
            rs = real(rs(abs(imag(rs))<epsilon & real(rs)>0 & real(rs)<arctimes(i)));
            % A local minimum or maximum of x(k+1) is found.
            if (~isempty(rs))
                for j = 1:length(rs)
                    xf = dynamics_onestep(x0(1:k+1),u,rs(j));
                    if (xf(k+1)<=M_min(k+2)-epsilon || xf(k+1)>=M_max(k+2)+epsilon)
                        flag_feasible = false;
                        return
                    end
                end
            end
        end
        x0 = dynamics_onestep(x0,u,arctimes(i));
        if (any(x0<=M_min(2:n+1)-epsilon | x0>=M_max(2:n+1)+epsilon))
            flag_feasible = false;
            return
        end
    end
    flag_feasible = true;
    return
end