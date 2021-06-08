  %Physical domain boundaries
  xa = 0.0;
  xb = 1.0;

  %Parameters
  smooth_opt   = ["none","weightavrg","gaussian"];
  avrg_opt     = ["arithmetic","harmonic","geometric"];
  interp_opt   = ["linear","pchip"];
  nnodes       = 201;
  time_start   = 0.0;
  velocity     = 1;
  diffusion    = 0.e-5;
  niter_mmpde  = 128;
  niter_coupl  = 2;
  eps_mmpde    = 1.e-3;
  smoothing    = smooth_opt(3);
  nsmooth      = 3;
  ksmooth      = 0.1;
  avrg         = avrg_opt(1); 
  numtimesteps = 0;
  timestep     = 0.0005;
  time 	       = time_start;
  eps_omega    = 1.e-2;
  interp_method= interp_opt(2);
  
  %Computational mesh is always uniform; ksi e [0,1]
  comp_mesh = linspace( 0, 1, nnodes)';
  
  %Initialize physical mesh as uniform
  physical_mesh = linspace( xa, xb, nnodes )';
  
  gamma_val = 0.5;



t0 = 0.0;
n_time_steps = 0;
tF = 0.2;
t = t0;

bl=1:3; br=nnodes-2:nnodes;

% Analytical solution
sol_exact_fun = @(x, t, vel) 0.0 + 1.0 * ( (x >= 0.25 + vel * t) & (x <= 0.5 + vel * t ) );

physical_grid = linspace( x0, x1, nnodes )';
uniform_grid = physical_grid;

velocity = 1.0;
cfl = 0.95;
dx_min = min(diff(physical_grid));

% initial condition evaluated on unifrom mesh
  %Boundary and initial conditions
  ua = 0.0;
  ub = 0.0;
solution = sol_exact_fun(physical_grid, 0.0, velocity) ; 

% compute the new time step
dt = cfl * dx_min / velocity;

omega = [];
physical_grid_uniform = [];

init_mesh = physical_grid;

while t < tF
    
    % Runge Kutta steps
%     
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
    
    %Solve nonlinear MMPDE to get final mesh
    [physical_mesh,~,u1,omega,l2res,itr_mmpde] = adapt_mesh(nnodes,xa,xb,comp_mesh,solution,solution,init_mesh,physical_mesh,niter_mmpde,eps_mmpde,eps_omega,avrg,smoothing,nsmooth,ksmooth,interp_method);

% % % % 
    solution = interp1(init_mesh,solution,physical_mesh,'pchip','extrap');
% 
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









