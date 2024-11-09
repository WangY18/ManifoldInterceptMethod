% Calculate the state vector at the end of the arc.
function xf = dynamics_onestep(x0, u, dt)
    n = length(x0);
    xf = zeros([n,1]);
    T = zeros([n+1,1]); % T(j+1) = dt^j/(j!)
    T(1) = 1;
    for j = 1:n
        T(j+1) = T(j)*dt/j;
    end
    for j = 1:n
        xf(j) = sum([x0(j:-1:1);u].*T(1:j+1));
    end
end