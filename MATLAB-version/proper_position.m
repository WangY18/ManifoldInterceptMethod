% Calculate the proper position of the trajectory.
function propos = proper_position(x0,xf,M_max,M_min,epsilon)
    n = length(x0);
    [orders,signs,tangents,arctimes] = plan_nth_order(x0(1:n-1),xf(1:n-1),M_max(1:n),M_min(1:n),false,0,epsilon);
    xs = end_points(x0,M_max(1),M_min(1),orders,signs,tangents,arctimes);
    propos = xs(n,end) - x0(n);
end