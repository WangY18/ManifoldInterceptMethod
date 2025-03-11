% Plan a 3rd-order trajectory from x0 to xf under constraint.
% First Step: 00, 010
% Second Step: s00, s010
% If direction is 1/-1, then directly consider the case where x0 is lower/higher than xf.
function [orders,signs,tangents,arctimes] = plan_3rd_order(x0,xf,M_max,M_min,flag_consider_position,direction,epsilon)
    if (direction==0)
        [orders_try,signs_try,tangents_try,arctimes_try] = plan_2nd_order(x0(1:2),xf(1:2),M_max(1:3),M_min(1:3),false,0,epsilon);
        if (isempty(orders_try))
            orders = [];
            signs = [];
            tangents = [];
            arctimes = [];
            return
        end
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
        if (xs_try(3,end)<xf(3))
            direction = 1;
        else
            direction = -1;
        end
    end

    if (direction<0) % x0 is higher than xf. It must be feasible.
        [orders,signs,tangents,arctimes] = plan_3rd_order(-x0,-xf,-M_min,-M_max,flag_consider_position,1,epsilon);
        signs = ~signs;
        return
    end

    x0 = x0(:);
    % x0 is lower than xf
    % Step 1: Try to enter the manifold.
    if (~isinf(M_max(3)))
        [orders,signs,tangents,arctimes] = plan_2nd_order(x0(1:2),[0;M_max(3)],M_max(1:3),M_min(1:3),false,1,epsilon);
        [orders2,signs2,tangents2,arctimes2] = plan_2nd_order([0;M_max(3)],xf(1:2),M_max(1:3),M_min(1:3),false,-1,epsilon);
        xs1 = end_points(x0,M_max(1),M_min(1),orders,signs,tangents,arctimes);
        xs2 = start_points(xf,M_max(1),M_min(1),orders2,signs2,tangents2,arctimes2);
        if (xs1(3,end)<=xs2(3,end))
            if (xs1(3,end)<xs2(3,end))
                orders = [orders,2];
                signs = [signs,true];
                tangents = [tangents,0];
                arctimes = [arctimes,(xs2(3,end)-xs1(3,end))/M_max(3)];
            end
            orders = [orders,orders2];
            signs = [signs,signs2];
            tangents = [tangents,tangents2];
            arctimes = [arctimes,arctimes2];
            return
        end
        xs1(:,end:-1:end-1) = [];
        orders(end) = [];
        signs(end) = [];
        tangents(end) = [];
        arctimes(end) = [];
        arctimes(end) = 0;
    elseif (~isinf(M_max(2)))
        orders = [0,1];
        signs = [true,true];
        tangents = [0,0];
        arctimes = [(M_max(2)-x0(1))/M_max(1),0];
        xs1 = [x0,dynamics_onestep(x0,M_max(1),arctimes(1))];
    else
        orders = 0;
        signs = true;
        tangents = 0;
        arctimes = 0;
        xs1 = x0;
    end
    % Step 2: Find where the trajectory enters the manifold.
    flag_succeed_noposition = false;
    propos0 = proper_position(x0,xf,M_max,M_min,epsilon);
    for i = size(xs1,2):-1:1
        propos = proper_position(xs1(:,i),xf,M_max,M_min,epsilon);
        if ((x0(3)+propos0-xf(3))*(xs1(3,i)+propos-xf(3))<=0)
            break;
        end
    end
    xs1(:,i+1:end) = [];
    orders(:,i+1:end) = [];
    signs(:,i+1:end) = [];
    tangents(:,i+1:end) = [];
    arctimes(:,i+1:end) = [];
    % xs1(:,end) enters the manifold along the arc.
    % Step 3: Create the whole trajectory. Consider s00, s010
    % Step 3.1: s00
    orders_try = [orders(end),0,0];
    signs_try = [signs(end),false,true];
    [T1,T2,T3] = solution_3arc_3rd_order(xs1(1:3,end),xf,orders_try,signs_try,M_max,M_min,epsilon);
    if (~isempty(T1))
        arctimes_try = [T1(1)-T2(1),T2(1)-T3(1),T3(1)];
        if (feasible(xs1(1:2,end),M_max(1:3),M_min(1:3),orders_try,signs_try,[0,0,0],arctimes_try,epsilon))
            flag_succeed_noposition = true;
            orders = [orders,orders_try(2:end)];
            signs = [signs,signs_try(2:end)];
            tangents = [tangents,0,0];
            arctimes = [arctimes(1:end-1),arctimes_try];
        end
    end
    if (~flag_succeed_noposition && ~isinf(M_min(2)))
        % Step 3.2: s010
        [order4,sign4,tangent4,arctime4] = plan_1st_order(M_min(2),xf(1),M_max(1),M_min(1));
        xsf = start_points(xf,M_max(1),M_min(1),order4,sign4,tangent4,arctime4);
        orders_try = [orders(end),0,1];
        signs_try = [signs(end),false,false];
        [T1,T2,T3] = solution_3arc_3rd_order(xs1(1:3,end),xsf(:,end),orders_try,signs_try,M_max,M_min,epsilon);
        if (~isempty(T1))
            orders_try = [orders_try,order4];
            signs_try = [signs_try,sign4];
            arctimes_try = [T1(1)-T2(1),T2(1)-T3(1),T3(1),arctime4];
            if (feasible(xs1(1:2,end),M_max(1:3),M_min(1:3),orders_try,signs_try,[0,0,0,0],arctimes_try,epsilon))
                orders = [orders,orders_try(2:end)];
                signs = [signs,signs_try(2:end)];
                tangents = [tangents,0,0,0];
                arctimes = [arctimes(1:end-1),arctimes_try];
                flag_succeed_noposition = true;
            end
        end
    end
    if flag_succeed_noposition && (~flag_consider_position || feasible(x0,M_max,M_min,orders,signs,tangents,arctimes,epsilon))
        return
    end
    % Step 4: Consider tangent markers.
    % Only -3 can support tangent markers, i.e., (-3,2).
    % Step 4.1: like 00(-3,2)000
    % Step 4.1.1: 00(-3,2)*
    orders = [0,0];
    signs = [true,false];
    [T1,T2] = solution_2arc_tangent_3rd_order(x0,orders,signs,M_max,M_min,true,epsilon);
    if (~isempty(T1))
        arctimes = [T1(1)-T2(1),T2(1)];
        tangents = [0,0];
        if (feasible(x0(1:2),M_max(1:3),M_min(1:3),orders,signs,tangents,arctimes,epsilon))
            xs1 = end_points(x0,M_max(1),M_min(1),orders,signs,tangents,arctimes);
            [orders2,signs2,tangents2,arctimes2] = plan_3rd_order(xs1(:,end),xf,M_max,M_min,false,-1,epsilon);
            orders = [orders,3,orders2];
            signs = [signs,true,signs2];
            tangents = [tangents,2,tangents2];
            arctimes = [arctimes,0,arctimes2];
            return
        end
    end
    % Step 4.1.2: 010(-3,2)*
    if (~isinf(M_max(2)))
        orders = [0,1,0];
        signs = [true,true,false];
        arctimes = [(M_max(2)-x0(1))/M_max(1),0,0];
        tangents = [0,0,0];
        x1 = dynamics_onestep(x0,M_max(1),arctimes(1));
        [T1,T2] = solution_2arc_tangent_3rd_order(x1,orders(2:3),signs(2:3),M_max,M_min,true,epsilon);
        if (~isempty(T1))
            arctimes(2:3) = [T1(1)-T2(1),T2(1)];
            if (feasible(x1(1:2),M_max(1:3),M_min(1:3),orders(2:3),signs(2:3),tangents(2:3),arctimes(2:3),epsilon))
                xs1 = end_points(x1,M_max(1),M_min(1),orders(2:3),signs(2:3),tangents(2:3),arctimes(2:3));
                [orders2,signs2,tangents2,arctimes2] = plan_3rd_order(xs1(:,end),xf,M_max,M_min,false,-1,epsilon);
                orders = [orders,3,orders2];
                signs = [signs,true,signs2];
                tangents = [tangents,2,tangents2];
                arctimes = [arctimes,0,arctimes2];
                return
            end
        end
    end
    % Step 4.2: like 000(-3,2)00
    % Step 4.2.1: *(-3,2)00
    orders = [0,0];
    signs = [false,true];
    [T1,T2] = solution_2arc_tangent_3rd_order(xf,orders,signs,M_max,M_min,false,epsilon);
    if (~isempty(T1))
        arctimes = [T2(1)-T1(1),-T2(1)];
        tangents = [0,0];
        xs1 = start_points(xf,M_max(1),M_min(1),orders(end:-1:1),signs(end:-1:1),tangents(end:-1:1),arctimes(end:-1:1));
        if (feasible(xs1(1:2,end),M_max(1:3),M_min(1:3),orders(end:-1:1),signs(end:-1:1),tangents(end:-1:1),arctimes(end:-1:1),epsilon))
            [orders2,signs2,tangents2,arctimes2] = plan_3rd_order(x0,xs1(:,end),M_max,M_min,false,-1,epsilon);
            orders = [orders2,3,orders(end:-1:1)];
            signs = [signs2,true,signs(end:-1:1)];
            tangents = [tangents2,2,tangents(end:-1:1)];
            arctimes = [arctimes2,0,arctimes(end:-1:1)];
            return
        end
    end
    % Step 4.2.1: *(-3,2)010
    if (~isinf(M_max(2)))
        orders = [0,1,0];
        signs = [false,true,true];
        arctimes = [(M_max(2)-x0(1))/M_max(1),0,0];
        x1 = dynamics_onestep(xf,M_min(1),-arctimes(1));
        [T1,T2] = solution_2arc_tangent_3rd_order(x1,orders(2:3),signs(2:3),M_max,M_min,false,epsilon);
        if (~isempty(T1))
            arctimes(1:2) = [T2(1)-T1(1),-T2(1)];
            tangents = [0,0,0];
            xs1 = start_points(xf,M_max(1),M_min(1),orders(end:-1:1),signs(end:-1:1),tangents(end:-1:1),arctimes(end:-1:1));
            if (feasible(xs1(1:2,end),M_max(1:3),M_min(1:3),orders(end:-1:1),signs(end:-1:1),tangents(end:-1:1),arctimes(end:-1:1),epsilon))
                [orders2,signs2,tangents2,arctimes2] = plan_3rd_order(x0,xs1(:,end),M_max,M_min,false,-1,epsilon);
                orders = [orders2,3,orders(end:-1:1)];
                signs = [signs2,true,signs(end:-1:1)];
                tangents = [tangents2,2,tangents(end:-1:1)];
                arctimes = [arctimes2,0,arctimes(end:-1:1)];
                return
            end
        end
    end
    % Between switching surfaces
    [orders,signs,tangents,arctimes] = plan_3rd_order(x0,xf,M_max,M_min,flag_consider_position,-direction,epsilon);
end