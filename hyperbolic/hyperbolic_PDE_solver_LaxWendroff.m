function [sol, tF] = hyperbolic_PDE_solver_LaxWendroff(x, sol, t0, c)

    n = length(x);
      
    u_old = sol;
    u_new = zeros(size(sol));
    
    cfl = 0.95;
    
    dx = diff(x);
    dx_min = min(dx);
    
    % compute the new time step
    dt = cfl * dx_min / c;
    
    for i = 2 : n - 1

       x_e = ( x(i)+x(i+1) )/2.0;
       x_w = ( x(i-1)+x(i) )/2.0;

       u_e = ( u_old(i)+u_old(i+1) )/2.0 - dt/(2*(x(i+1)-x(i))) * c * ( u_old(i+1) - u_old(i) );
       u_w = ( u_old(i-1)+u_old(i) )/2.0 - dt/(2*(x(i)-x(i-1))) * c * ( u_old(i) - u_old(i-1) );

       u_new(i) = u_old(i) -  dt/(x_e-x_w) * c * ( u_e - u_w ); 

    end
    
    sol = u_new; 
    
    tF = t0 + dt;

return