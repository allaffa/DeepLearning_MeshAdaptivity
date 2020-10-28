function [value] = speed_sound(gamma, p, rho)

    a = sqrt(gamma * p./rho);
    value = max(a);

end