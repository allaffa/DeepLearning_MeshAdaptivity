function [x, u1_new, u2_new, u3_new, p_new, t] = Euler_equations_solver_conservative_WENO5(gamma_var, x, U1, U2, P, t0, tF, nghost)

    n = length(x);
    dx_min = min(diff(x));
%     bl=1:3; br=n-2:n;
  
    % initial condition for density
    u1_old = U1;
    
    % initial condition for velocity
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
               
    %% Lax-Wendroff 

    while(t<tF)
          
        % compute the new time step
        dt = time_step_compute(gamma_var, p_old, u1_old, u2_old./u1_old, x);
        
        %% Runge-Kutta step    
        q0 = sol; 
        q0=refletive_bc_apply(q0,nghost);
                
        % 1st stage
        res= FV_compWise_WENO5LF1d(gamma_var,lambda,dt,q0,dx_min,nghost);
        k1 = -res;
        q1 = q0 + 0.5*k1*dt;
        q1=refletive_bc_apply(q1,nghost);
%         q=q0-dt*res;
%         q(:,bl)=q0(:,bl); q(:,br)=q0(:,br); % Neumann BCs

        % 2nd Stage
        res=FV_compWise_WENO5LF1d(gamma_var,lambda,0.5*dt,q1,dx_min,nghost);
        k2 = -res;
        q2 = q0 + ( -k1 + 2.0*k2 ) * dt;
        q2=refletive_bc_apply(q2,nghost);
%         q = 0.75*q0+0.25*(q-dt*res);
%         q(:,bl)=q0(:,bl); q(:,br)=q0(:,br); % Neumann BCs

        % 3rd stage
        res=FV_compWise_WENO5LF1d(gamma_var,lambda,dt,q2,dx_min,nghost);
        k3 = -res;
        q = q0 + ( k1 + 4.0*k2 + k3 ) / 6.0 * dt;
        q=refletive_bc_apply(q,nghost);
%         q = (q0+2*(q-dt*res))/3;
%         q(:,bl)=q0(:,bl); q(:,br)=q0(:,br); % Neumann BCs

        sol = q;     
        
        u1_new = sol(1,:);
        u2_new = sol(2,:);
        u3_new = sol(3,:);
        p_new = (gamma_var - 1) * ( u3_new - 1/2 * u2_new .* u2_new ./ u1_new );
        
        u1_old = u1_new;
        u2_old = u2_new;
        u3_old = u3_new;
        p_old = p_new;
        
        a=sqrt(gamma_var*p_old./u1_old); 

        % Compute dt for next time step
        lambda=max(abs(u2_old./u1_old)+a);         
        
        t = t + dt;
        
    end  
   
end
   
