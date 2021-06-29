function [x, u1_new, u2_new, u3_new, p_new, t, time_history, grid_history] = Euler_equations_solver_conservative_WENO5_nonuniform(gamma_var, x, U1, U2, P, t0, tF, nghost, verbose, kernel_width, adaptive_mesh, python_module)

    n = length(x);

    time_history = [];
    grid_history = [];
    
    dx = zeros(size(x));
    dx(1) = x(2)-x(1);
    dx(end) = x(end) - x(end-1);

    for i =2:length(x)-1
        dx(i) = 0.5 * ((x(i)-x(i-1))+(x(i+1)-x(i)));
    end

    bl=1:3; br=n-2:n;
  
    % initial condition for density
    u1_old = U1;
    
    % initial condition for density * velocity
    u2_old = U2;
    
    %initial condition for pressure
    p_old = P; 
    
    % initial mesh
    x_old = x;
    
    % total energy is derived from the pressure using e.o.s.
    u3_old = p_old./(gamma_var-1) + 1/2 * u2_old .* u2_old ./ u1_old;
    
    sol = [U1'; U2'; u3_old'];    
    
    u1_new = zeros(size(u1_old));
    u2_new = zeros(size(u2_old));
    u3_new = zeros(size(u3_old));
    p_new = zeros(size(p_old));
    
    t = t0;
    
    a=sqrt(gamma_var*p_old./u1_old); 
    
    % Compute dt for next time step
    lambda=max(abs(u2_old./u1_old)+a); 
    
    count = 1;
               
    %% Lax-Wendroff 

    while(t<tF)
          
        uniform_mesh = linspace( x(1), x(end), n)';

        %Parameters
        smooth_opt   = ["none","weightavrg","gaussian"];
        avrg_opt     = ["arithmetic","harmonic","geometric"];
        interp_opt   = ["linear","pchip","spline"];
        niter        = 50;
        eps_nonlin   = 1.e-5;
        smoothing    = smooth_opt(1);
        nsmooth      = 10;
        ksmooth      = 0.1;
        avrg         = avrg_opt(1); 
        interp_method= interp_opt(2);

        if adaptive_mesh == "standard" && mod(count, 1)==0   

           eps_omega    = 1.e-8;
           [x_new,~,~,~,~,~] = adapt_mesh(n-2*nghost,x(nghost+1),x(end-nghost),uniform_mesh(nghost+1:end-nghost),u2_old(nghost+1:end-nghost),u2_old(nghost+1:end-nghost),x_old(nghost+1:end-nghost),x_old(nghost+1:end-nghost),niter,eps_nonlin,eps_omega,avrg,smoothing,nsmooth,ksmooth,interp_method);

           x_new = [x(1:nghost); x_new; x(end-nghost+1:end)];
           time_history = [time_history; t];
           grid_history = [grid_history; x_new'];

           u1_old2 = interp1(x_old,u1_old,x_new, 'pchip', 'extrap');
           u2_old2 = interp1(x_old,u2_old,x_new, 'pchip', 'extrap');
           u3_old2 = interp1(x_old,u3_old,x_new, 'pchip', 'extrap');
           p_old2 = interp1(x_old,p_old,x_new, 'pchip', 'extrap');   

           %Initial condition for the next time step
           u1_old = u1_old2;
           u2_old = u2_old2;
           u3_old = u3_old2;  
           p_old = p_old2;
           x = x_new;
           x_old = x_new;  
       
           dx = zeros(size(x));
           dx(1) = x(2)-x(1);
           dx(end) = x(end) - x(end-1);

           for i =2:length(x)-1
               dx(i) = 0.5 * ((x(i)-x(i-1))+(x(i+1)-x(i)));
           end

           sol(1,:) = u1_old;
           sol(2,:) = u2_old;
           sol(3,:) = u3_old;

        elseif adaptive_mesh == "DL" && mod(count, 1)==0
            
           eps_omega    = 1.e-6;
           grad_u = gradient(u1_old(nghost+1:end-nghost),x_old(nghost+1:end-nghost));
           monitor_function   = monitor_fun(grad_u,eps_omega);
           monitor_function = interp1(x_old(nghost+1:end-nghost),monitor_function,uniform_mesh(nghost+1:end-nghost),interp_method,'extrap');
           AI_mesh = double(py.array.array('d',py.numpy.nditer(python_module.evaluate_model_load(monitor_function))));
           AI_mesh = cumsum(AI_mesh);
           x_new = x(nghost+1) + x_old(end-nghost)*(AI_mesh - AI_mesh(1))/(AI_mesh(end) - AI_mesh(1));
           x_new = x_new';

           x_new = [x(1:nghost); x_new; x(end-nghost+1:end)];
           time_history = [time_history; t];
           grid_history = [grid_history; x_new'];

           u1_old2 = interp1(x_old,u1_old,x_new, 'spline', 'extrap');
           u2_old2 = interp1(x_old,u2_old,x_new, 'spline', 'extrap');
           u3_old2 = interp1(x_old,u3_old,x_new, 'spline', 'extrap');
           p_old2 = interp1(x_old,p_old,x_new, 'spline', 'extrap');   

           %Initial condition for the next time step
           u1_old = u1_old2;
           u2_old = u2_old2;
           u3_old = u3_old2;  
           p_old = p_old2;
           x = x_new;
           x_old = x_new;     
      
           dx = zeros(size(x));
           dx(1) = x(2)-x(1);
           dx(end) = x(end) - x(end-1);

           for i =2:length(x)-1
               dx(i) = 0.5 * ((x(i)-x(i-1))+(x(i+1)-x(i)));
           end

           sol(1,:) = u1_old;
           sol(2,:) = u2_old;
           sol(3,:) = u3_old;

        end

        % LUKA MAY ASK ABOUT THIS
        u1_old2 = interp1(x_old,u1_old,x_old, 'pchip', 'extrap');
        u2_old2 = interp1(x_old,u2_old,x_old, 'pchip', 'extrap');
        u3_old2 = interp1(x_old,u3_old,x_old, 'pchip', 'extrap');
        p_old2 = interp1(x_old,p_old,x_old, 'pchip', 'extrap');   

        %Initial condition for the next time step
        u1_old = u1_old2;
        u2_old = u2_old2;
        u3_old = u3_old2;  
        p_old = p_old2;   

        % compute the new time step
        dt = time_step_compute(gamma_var, p_old, u1_old, u2_old./u1_old, x);
        
        if(verbose)
           display(strcat('Time(s): ', num2str(t)));
        end        
        
        %% Runge-Kutta step    
        q0 = sol; 
                
        % 1st stage
        res= FV_compWise_WENO5LF1d_nonuniform(gamma_var,lambda,dt,q0,dx);
        k1 = -res;
        q1 = q0 + 0.5*k1*dt;
        q1=refletive_bc_apply(q1,nghost);

        % 2nd Stage
        res=FV_compWise_WENO5LF1d_nonuniform(gamma_var,lambda,dt,q1,dx);
        k2 = -res;
        q2 = q0 + ( -k1 + 2.0*k2 ) * dt;
        q2=refletive_bc_apply(q2,nghost);

        % 3rd stage
        res=FV_compWise_WENO5LF1d_nonuniform(gamma_var,lambda,dt,q2,dx);
        k3 = -res;
        q = q0 + ( k1 + 4.0*k2 + k3 ) / 6.0 * dt;
        q=refletive_bc_apply(q,nghost);

        sol = q;     
        
        u1_new = sol(1,:);
        u2_new = sol(2,:);
        u3_new = sol(3,:);
        p_new = (gamma_var - 1) * ( u3_new - 1/2 * u2_new .* u2_new ./ u1_new );
        
        u1_old = u1_new';
        u2_old = u2_new';
        u3_old = u3_new';
        p_old = p_new';
        
        a=sqrt(gamma_var*p_old./u1_old); 

        % Compute dt for next time step
        lambda=max(abs(u2_old./u1_old)+a);         
        
        t = t + dt;
        count = count+1;
        
    end  
   
end
   
