x0 = 2;
x1 = 4;
nnodes = 500;

t0 = 0.0;
n_time_steps = 50;
dt = 0.005;
t = t0;

kernel_width = nnodes * 0.000001;

% Analytical solution
sol_exact_fun = @(x, t, vel) 0.0 + 1.0 * ( (x >= 2.5 + vel * t) & (x <= 3 + vel * t ) );

velocity = 1.0;

% time integration (Explicit Eulfer: alpha=0.0; Implicit Eulfer: alpha=1.0; Crank-Nicholson: alpha = 0.5)
alpha = 0.5;

physical_grid = linspace( x0, x1, nnodes )';

% initial condition evaluated on unifrom mesh
solution = sol_exact_fun(physical_grid, 0.0, velocity) ; 

% initial condition evaluated on adapted mesh
% solution = sol_exact_fun(physical_grid, 0.0, velocity) ; 

for time_step = 1:n_time_steps 
    
    new_physical_grid = adaptive_mesh_moving1D_optimization(nnodes, physical_grid, solution, kernel_width);

    % interpolate solution from previous physical grid to new physical grid
    solution(2:end-1) = interp1(physical_grid, solution, new_physical_grid(2:end-1));    
    
    physical_grid = new_physical_grid;
    
    solution = hyperbolic_PDE_solver(physical_grid, solution, dt, velocity, alpha);
    
    t = t + dt;
        
end

%figure()
plot(physical_grid, sol_exact_fun(physical_grid, t, velocity), 'linewidth', 2);
set(gca, 'fontsize', 40);
hold on
plot(physical_grid, solution, '-', 'linewidth', 3);
%     legend('Analytical solution', 'Numerical solution');
title(strcat('Time: ', num2str(t)))

% figure()
% plot(linspace( x0, x1, nnodes ), physical_grid, '*', 'LineWidth', 5)



