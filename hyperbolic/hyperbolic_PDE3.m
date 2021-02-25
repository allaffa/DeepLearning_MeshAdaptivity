clear

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

% close all
clc

x0 = -10;
x1 = 10;
nnodes = 201;

t0 = 0.0;
n_time_steps = 100;
dt = 0.005;
t = t0;

kernel_width = nnodes * 0.000001;

% Analytical solution
sol_exact_fun = @(x, t, vel) exp(-((x + 0.5 - vel * t).^2)/0.01);

velocity = 1.0;

% time integration (Explicit Eulfer: alpha=0.0; Implicit Eulfer: alpha=1.0; Crank-Nicholson: alpha = 0.5)
alpha = 0.5;

physical_grid = linspace( x0, x1, nnodes )';
uniform_grid = physical_grid;

% initial condition evaluated on unifrom mesh
solution = sol_exact_fun(physical_grid, 0.0, velocity) ; 

solution_uniform = solution;

for time_step = 1:n_time_steps 
    
    new_physical_grid = adaptive_mesh_moving1D_optimization(physical_grid, solution, kernel_width);
    
    if time_step == 1
        solution = sol_exact_fun(new_physical_grid, 0.0, velocity) ; 
    else
        % interpolate solution from previous physical grid to new physical grid
        solution(2:end-1) = interp1(physical_grid, solution, new_physical_grid(2:end-1));
%         solution(2:end-1) = pchip(physical_grid, solution, new_physical_grid(2:end-1));        
    end
    
    physical_grid = new_physical_grid;
    
%     solution = hyperbolic_PDE_solver(physical_grid, solution, dt, velocity, alpha);  
%     
%     solution_uniform = hyperbolic_PDE_solver(uniform_grid, solution_uniform, dt, velocity, alpha);
    
    t = t + dt;
        
end

figure()
plot(linspace( x0, x1, nnodes ), physical_grid, '-', 'LineWidth', 10)
hold on 
axis equal
xlim([x0 x1])
ylim([x0 x1])
xlabel('Uniform mesh coordinates')
ylabel('Non-uniform mesh coordinates')
legend('standard adaptivity', 'deep learning')
set(gca, 'fontsize', 45);

figure()
plot(physical_grid, sol_exact_fun(physical_grid, t, velocity), 'linewidth', 10);
xlim([x0 x1])
ylim([-0.2 1.2])
set(gca, 'fontsize', 45);
xlabel('Coordinates in the physical domain')
ylabel('Solution values')
hold on
plot(uniform_grid, solution_uniform, 'o-', 'linewidth', 10);
plot(physical_grid, solution, 'o-', 'linewidth', 10);
    legend('Analytical solution', 'Numerical - uniform mesh', 'Numerical - adaptive mesh');
title(strcat('Time: ', num2str(t)))




