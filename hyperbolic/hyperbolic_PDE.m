close all
clear

%%

pyExec	= '/opt/anaconda3/bin/python3';
% [ver, exec, loaded]	= pyversion(pyExec);

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

x0 = 0;
x1 = 1;
nnodes = 201;

t0 = 0.00;
tF = 0.05;
t = t0;
alpha = 1.0;

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
    
    % compute the new time step
    dt = cfl * dx_min / velocity;
    
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


%         [new_physical_grid,~,~,monr_function,~,~] = adapt_mesh(nnodes,x0,x1,uniform_grid, solution,solution,physical_grid,physical_grid,niter,eps_nonlin,eps_omega,avrg,smoothing,nsmooth,ksmooth,interp_method);
% 
           grad_u = gradient(solution,physical_grid);
           monitor_function   = monitor_fun(grad_u,eps_omega);
           monitor_function = interp1(physical_grid,monitor_function,uniform_grid,interp_method,'extrap');
           AI_mesh = double(py.array.array('d',py.numpy.nditer(python_module.evaluate_model_load(monitor_function))));
           AI_mesh = cumsum(AI_mesh);
           new_physical_grid = x0 + (AI_mesh - AI_mesh(1))/(AI_mesh(end) - AI_mesh(1));
           new_physical_grid = new_physical_grid';

        solution(2:end-1) = interp1(physical_grid, solution, new_physical_grid(2:end-1),'pchip','extrap');
        
        physical_grid = new_physical_grid;
        
    end
    
    solution = hyperbolic_PDE_solver(physical_grid, solution, dt, velocity, alpha);
    
    t = t + dt;
        
end

t = t-dt;


figure()
plot(linspace( x0, x1, nnodes ), physical_grid, '-', 'LineWidth', 10)
hold on 
axis equal
xlim([x0 x1])
ylim([x0 x1])
xlabel('Coordinates of uniform mesh')
ylabel('Coordinates of adapted mesh')
legend('standard adaptivity', 'deep learning')
set(gca, 'fontsize', 45);


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









