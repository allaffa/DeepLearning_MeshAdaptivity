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

t0 = 0.0;
tF = 0.4;
t = t0;

%Parameters
smooth_opt   = ["none","weightavrg","gaussian"];
avrg_opt     = ["arithmetic","harmonic","geometric"];
interp_opt   = ["linear","pchip","spline"];
time_start   = 0.0;
velocity     = 1.0;
diffusion    = 0.e-5;
niter_mmpde  = 50;
niter_coupl  = 2;
eps_mmpde    = 1.e-3;
smoothing    = smooth_opt(1);
nsmooth      = 10;
ksmooth      = 0.12;
avrg         = avrg_opt(1); 
numtimesteps = 800;
timestep     = 0.0005;
time 	       = time_start;
eps_omega    = 0.5*1e-3;
interp_method= interp_opt(2);

adaptivity = true;

bl=1:3; br=nnodes-2:nnodes;

% Analytical solution
sol_exact_fun = @(x, t, vel) 0.0 + 1.0 * ( ( (x >= 0.30 + vel * t) & (x <= 0.50 + vel * t ) ) ) ;

physical_grid = linspace( x0, x1, nnodes )';
uniform_grid = physical_grid;

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

%         [new_physical_grid,~,~,monr_function,~,~] = adapt_mesh(nnodes,x0,x1,uniform_grid, solution,solution,physical_grid,physical_grid,niter_mmpde,eps_mmpde,eps_omega,avrg,smoothing,nsmooth,ksmooth,interp_method);
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
plot(max(0,physical_grid), sol_exact_fun(physical_grid, t, velocity), 'linewidth', 6);
set(gca, 'fontsize', 72);
hold on
ylim([-0.2, 1.2])
plot(max(0,physical_grid), solution, '-', 'linewidth', 6);
    legend('Exact solution', 'Numerical solution');
xlabel('Coordinates in the physical domain')
ylabel('Solution values')
title(strcat('Time: ', num2str(t)))

figure()
plot(diff(physical_grid), '-o', 'Linewidth', 4)
title('Mesh spacing')
xlabel('Node index')
ylabel('\Delta x')
set(gca, 'fontsize', 45)









