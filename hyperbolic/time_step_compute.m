function [value] = time_step_compute(gamma, p, rho, velocity, x)

    a = speed_sound(gamma, p, rho);
    cfl = 0.55;
    dx = diff(x);
    dx_min = min(dx);
    u_max = max(abs(velocity));
    
    value = cfl .* dx_min/(u_max + a);

end