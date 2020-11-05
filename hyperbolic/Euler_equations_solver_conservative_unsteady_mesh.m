function [x, u1_new, u2_new, u3_new, p_new, t, monitor_function] = Euler_equations_solver_conservative_unsteady_mesh(gamma, x, U1, U2, P, t0, tF, upwind, verbose, solver, kernel_width, adaptive_mesh)

    n = length(x);

    % Explicit time stepping (Forward Euler)
    
    % Convention: u1 = density
    %             u2 = velocity
    %             u3 = energy
    %             p = pressure
    
    % Hyperbolic equation: du/dt + f(u) = 0

    % nonlinear term for density
    F1 = @(rho, rho_u, E, p) rho_u; 
    
    %nonlinear term for velocity
    %F2 = @(rho, u, E, p) u .* u + p./rho;
    F2 = @(rho, rho_u, E, p) rho_u .* rho_u ./ rho + p;
    
    %nonlinear term for energy
    F3 = @(rho, rho_u, E, p) rho_u./rho .* (E+p); 
    
    % initial condition for density
    u1_old = U1;
    
    % initial condition for velocity
    u2_old = U2;
    
    %initial condition for pressure
    p_old = P; 
    
    % initial mesh
    x_old = x;
    
    % total energy is derived from the pressure using e.o.s.
    u3_old = p_old./(gamma-1) + 1/2 * u2_old .* u2_old ./ u1_old;
    
    u1_new = zeros(size(u1_old));
    u2_new = zeros(size(u2_old));
    u3_new = zeros(size(u3_old));
    p_new = zeros(size(p_old));
    
    t = t0;
    
           
   %% Lax-Wendroff 
   
   while(t<tF)

       % compute the new time step
       dt = time_step_compute(gamma, p_old, u1_old, u2_old./u1_old, x);       
       
       if(adaptive_mesh)
           
           [x_new, ~, monitor_function] = unsteady_adaptive_mesh_moving1D_optimization( x_old, u1_old, kernel_width, 1e-7 );
           u1_old2 = interp1(x_old,u1_old,x_new, 'linear', 'extrap');
           u2_old2 = interp1(x_old,u2_old,x_new, 'linear', 'extrap');
           u3_old2 = interp1(x_old,u3_old,x_new, 'linear', 'extrap');
           p_old2 = interp1(x_old,p_old,x_new, 'linear', 'extrap');   

           %Initial condition for the next time step
           u1_old = u1_old2;
           u2_old = u2_old2;
           u3_old = u3_old2;  
           p_old = p_old2;
           x = x_new;
           x_old = x_new;
           
           % Re-compute the new time step
           dt = time_step_compute(gamma, p_old, u1_old, u2_old./u1_old, x);           
           
       end       

       if(verbose)
           display(strcat('Time(s): ', num2str(t)));
       end
       % For each time step, update each element in the grid
       for i = 2 : n - 1

           x_e = ( x(i)+x(i+1) )/2.0;
           x_w = ( x(i-1)+x(i) )/2.0;

           F1_E = F1(u1_old(i+1), u2_old(i+1), u3_old(i+1), p_old(i+1));
           F1_P = F1(u1_old(i), u2_old(i), u3_old(i), p_old(i));
           F1_W = F1(u1_old(i-1), u2_old(i-1), u3_old(i-1), p_old(i-1));

           F2_E = F2(u1_old(i+1), u2_old(i+1), u3_old(i+1), p_old(i+1));
           F2_P = F2(u1_old(i), u2_old(i), u3_old(i), p_old(i));
           F2_W = F2(u1_old(i-1), u2_old(i-1), u3_old(i-1), p_old(i-1));

           F3_E = F3(u1_old(i+1), u2_old(i+1), u3_old(i+1), p_old(i+1));
           F3_P = F3(u1_old(i), u2_old(i), u3_old(i), p_old(i));
           F3_W = F3(u1_old(i-1), u2_old(i-1), u3_old(i-1), p_old(i-1));

           % DOUBLE CHECK THAT THIS MAKES SENSE WHEN THE MESH IS NOT
           % UNIFORM
           u1_e = ( u1_old(i)+u1_old(i+1) )/2.0 - dt/(2*(x(i+1)-x(i))) * ( F1_E - F1_P );
           u1_w = ( u1_old(i-1)+u1_old(i) )/2.0 - dt/(2*(x(i)-x(i-1))) * ( F1_P - F1_W );

           u2_e = ( u2_old(i)+u2_old(i+1) )/2.0 - dt/(2*(x(i+1)-x(i))) * ( F2_E - F2_P );
           u2_w = ( u2_old(i-1)+u2_old(i) )/2.0 - dt/(2*(x(i)-x(i-1))) * ( F2_P - F2_W );

           u3_e = ( u3_old(i)+u3_old(i+1) )/2.0 - dt/(2*(x(i+1)-x(i))) * ( F3_E - F3_P );
           u3_w = ( u3_old(i-1)+u3_old(i) )/2.0 - dt/(2*(x(i)-x(i-1))) * ( F3_P - F3_W );

           p_e = (gamma - 1) * ( u3_e - 1/2 * u2_e .* u2_e ./ u1_e );
           p_w = (gamma - 1) * ( u3_w - 1/2 * u2_w .* u2_w ./ u1_w );

           F1_e = F1(u1_e, u2_e, u3_e, p_e);
           F1_w = F1(u1_w, u2_w, u3_w, p_w);

           F2_e = F2(u1_e, u2_e, u3_e, p_e);
           F2_w = F2(u1_w, u2_w, u3_w, p_w);

           F3_e = F3(u1_e, u2_e, u3_e, p_e);
           F3_w = F3(u1_w, u2_w, u3_w, p_w);

           u1_new(i) = u1_old(i) - dt/(x_e-x_w) * ( F1_e - F1_w );
           u2_new(i) = u2_old(i) - dt/(x_e-x_w) * ( F2_e - F2_w );
           u3_new(i) = u3_old(i) - dt/(x_e-x_w) * ( F3_e - F3_w );

           p_new(i) = (gamma - 1) * ( u3_new(i) - 1/2 * u2_new(i) .* u2_new(i) ./ u1_new(i) );

       end
        
       %Boundary values extrapolation
       % Not done now
       u1_new(1) = u1_new(2);
       u2_new(1) = u2_new(2);
       u3_new(1) = u3_new(2);
       p_new(1) = p_new(2);
       u1_new(size(u1_new,1)) = u1_new(size(u1_new,1)-1);
       u2_new(size(u2_new,1)) = u2_new(size(u2_new,1)-1);
       u3_new(size(u3_new,1)) = u3_new(size(u3_new,1)-1);
       p_new(size(p_new,1)) = p_new(size(p_new,1)-1); 

       %Initial condition for the next time step
       u1_old = u1_new;
       u2_old = u2_new;
       u3_old = u3_new;  
       p_old = p_new; 

       t = t + dt;

   end   
   
   
end