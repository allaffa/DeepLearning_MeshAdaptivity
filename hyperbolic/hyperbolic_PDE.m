pyExec	= '/opt/anaconda3/bin/python3';

python_path = py.sys.path;

% Add folders to python system path to loead version 3.7 of Python.
if count(python_path, pyExec) == 0
    insert(py.sys.path, int64(0), pyExec);
end

% Verify that the Python 3.7 version from Anaconda is successfully loaded
[ver2, exec2, loaded2]	= pyversion;
assert(loaded2==1);

python_module = py.importlib.import_module('model_load_script');

%%

close all
clc

x0 = 0;
x1 = 1;
nnodes = 501;

t0 = 0.0;
n_time_steps = 60;
dt = 0.000005;
tF = 0.2;
t = t0;

kernel_width = nnodes * 0.00000001;

% Analytical solution
sol_exact_fun = @(x, t, vel) 0.0 + 1.0 * ( (x >= 0.25 + vel * t) & (x <= 0.5 + vel * t ) );

velocity = 1.0;

% time integration (Explicit Eulfer: alpha=0.0; Implicit Eulfer: alpha=1.0; Crank-Nicholson: alpha = 0.5)
alpha = 0.0;

physical_grid = linspace( x0, x1, nnodes )';
uniform_grid = physical_grid;

% initial condition evaluated on unifrom mesh
solution = sol_exact_fun(physical_grid, 0.0, velocity) ; 

% initial condition evaluated on adapted mesh
% solution = sol_exact_fun(physical_grid, 0.0, velocity) ; 

omega = [];
physical_grid_uniform = [];

while t < tF
    
%     [new_physical_grid, physical_grid_uniform, omega] = adaptive_mesh_moving1D_optimization(physical_grid, solution, kernel_width);

    % interpolate solution from previous physical grid to new physical grid
%     solution(2:end-1) = interp1(physical_grid, solution, new_physical_grid(2:end-1));    
    
%     physical_grid = new_physical_grid;
    
    solution = hyperbolic_PDE_solver(physical_grid, solution, dt, velocity, alpha);
    
    t = t + dt;
        
end


% figure()
% plot(linspace( x0, x1, nnodes ), physical_grid, '*', 'LineWidth', 10)
% hold on 
% axis equal
% xlim([x0 x1])
% ylim([x0 x1])
% xlabel('Uniform mesh coordinates')
% ylabel('Non-uniform mesh coordinates')
% legend('standard adaptivity', 'deep learning')
% set(gca, 'fontsize', 45);


figure()
plot(max(0,physical_grid), sol_exact_fun(physical_grid, t, velocity), 'linewidth', 5);
set(gca, 'fontsize', 45);
ylim([-0.2, 1.2])
hold on
% plot(max(0,physical_grid), 0.5, 'o', 'Linewidth', 3)
% plot(max(0,physical_grid_uniform), omega, '-o', 'Linewidth', 1)
% %ylim([-0.2 1.2])
% xlabel('Coordinates in the physical domain')
% ylabel('Solution values')
plot(max(0,physical_grid), solution, '-', 'linewidth', 5);
    legend('Analytical solution', 'Numerical solution');
title(strcat('Time: ', num2str(t)))









