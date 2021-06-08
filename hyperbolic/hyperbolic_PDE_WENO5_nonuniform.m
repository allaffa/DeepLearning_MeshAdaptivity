% close all
% clc

x0 = 0;
x1 = 1;
nnodes = 501;

t0 = 0.0;
n_time_steps = 100;
tF = 0.2;
t = t0;

adaptivity = false;

bl=1:3; br=nnodes-2:nnodes;

kernel_width = nnodes * 0.00000001;

% Analytical solution
sol_exact_fun = @(x, t, vel) 0.0 + 1.0 * ( (x >= 0.25 + vel * t) & (x <= 0.5 + vel * t ) );

physical_grid = linspace( x0, x1, nnodes )';
uniform_grid = physical_grid;

velocity = 1.0;
cfl = 0.2;
dx_min = min(diff(physical_grid));

% initial condition evaluated on unifrom mesh
solution = sol_exact_fun(physical_grid, 0.0, velocity) ; 

% compute the new time step
dt = cfl * dx_min / velocity;

omega = [];
physical_grid_uniform = [];

while t < tF
    
    if (adaptivity)

        % Runge Kutta steps
    %     
          %Parameters
          smooth_opt   = ["none","weightavrg","gaussian"];
          avrg_opt     = ["arithmetic","harmonic","geometric"];
          interp_opt   = ["linear","pchip","spline"];
          niter        = 50;
          eps_nonlin   = 1.e-5;
          smoothing    = smooth_opt(2);
          nsmooth      = 10;
          ksmooth      = 0.2;
          avrg         = avrg_opt(1); 
          eps_omega    = 1.e-2;
          interp_method= interp_opt(2);


        [new_physical_grid,~,~,monr_function,~,~] = adapt_mesh(nnodes,x0,x1,uniform_grid, solution,solution,physical_grid,physical_grid,niter,eps_nonlin,eps_omega,avrg,smoothing,nsmooth,ksmooth,interp_method);
        
        solution(2:end-1) = interp1(physical_grid, solution, new_physical_grid(2:end-1),'pchip','extrap');
        
        physical_grid = new_physical_grid;
        
    end
    
    q0 = solution; 
       
    % 1st stage
    res= hyperbolic_PDE_solver_WENO5_nonuniform(physical_grid, q0, velocity, dt);
    q=q0-dt*res;
    q(bl)=q0(bl); q(br)=q0(br); % Neumann BCs
        
    % 2nd Stage
    res=hyperbolic_PDE_solver_WENO5_nonuniform(physical_grid, q, velocity, dt);  
    q = 0.75*q0+0.25*(q-dt*res);
    q(bl)=q0(bl); q(br)=q0(br); % Neumann BCs
    
    % 3rd stage
    res=hyperbolic_PDE_solver_WENO5_nonuniform(physical_grid, q, velocity, dt);   
    q = (q0+2*(q-dt*res))/3;
    q(bl)=q0(bl); q(br)=q0(br); % Neumann BCs

    solution = q;
    
    t = t + dt;
end


figure()
plot(linspace( x0, x1, nnodes ), physical_grid, '*', 'LineWidth', 10)
hold on 
axis equal
xlim([x0 x1])
ylim([x0 x1])
xlabel('Uniform mesh coordinates')
ylabel('Non-uniform mesh coordinates')
legend('standard adaptivity', 'deep learning')
set(gca, 'fontsize', 45);


figure()
plot(max(0,physical_grid), sol_exact_fun(physical_grid, t, velocity), 'linewidth', 4);
set(gca, 'fontsize', 45);
ylim([-0.2, 1.2])
hold on
% plot(max(0,physical_grid), 0.5, 'o', 'Linewidth', 3)
% plot(max(0,physical_grid_uniform), omega, '-o', 'Linewidth', 1)
% %ylim([-0.2 1.2])
% xlabel('Coordinates in the physical domain')
% ylabel('Solution values')
plot(max(0,physical_grid), solution, '-', 'linewidth', 4);
    legend('Analytical solution', 'Numerical solution');
title(strcat('Time: ', num2str(t)))









