function [xs,ts,Tnext] = interpolate_MIM(x0,orders,signs,tangents,arctimes,M0_max,M0_min,Ts,T0,flag_end)
    if nargin<9
        T0 = 0;
    end
    if nargin<10
        flag_end = true;
    end
    toltime = sum(arctimes(tangents==0));
    x0 = x0(:);
    T0_ = T0;
    n = length(x0);
    xs = zeros([n+1,floor((toltime-T0)/Ts)+1]);
    idx = 1;
    for i = 1:length(orders)
        if (tangents(i)~=0)
            continue
        end
        arctime = arctimes(i);
        u = cal_u(orders(i),signs(i),M0_max,M0_min);
        while (T0<=arctime)
            xs(1,idx) = u;
            xs(2:n+1,idx) = dynamics_onestep(x0,u,T0);
            idx = idx + 1;
            T0 = T0 + Ts;
        end
        T0 = T0 - arctime;
        x0 = dynamics_onestep(x0,u,arctime);
    end
    ts = T0_ + Ts*(0:length(xs)-1);
    if (flag_end)
        xs = [xs,[u;x0]];
        ts = [ts,toltime];
    end
    Tnext = Ts - T0;
end