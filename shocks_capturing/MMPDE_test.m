clear all
%close all
%clc

  %Problem definiton
  %exact_sol = @(x, t, v) 0.0 + 1.0 * ( ( (x >= 0.25 + v * t) & (x <= 0.45 + v * t ) ) ) + 1.0 * ( ( (x >= 0.65 + v * t) & (x <= 0.85 + v * t ) ) );
  exact_sol = @(x, t, v) 0.0 + 1.0 * ( ( (x >= 0.30 + v * t) & (x <= 0.50 + v * t ) ) ) ;
  
  %Physical domain boundaries
  xa = 0.2;
  xb = 0.9;

  %Parameters
  smooth_opt   = ["none","weightavrg","gaussian"];
  avrg_opt     = ["arithmetic","harmonic","geometric"];
  interp_opt   = ["linear","pchip","spline"];
  nnodes       = 81;
  time_start   = 0.0;
  velocity     = 1.e-2;
  diffusion    = 0.e-5;
  niter_mmpde  = 64;
  niter_coupl  = 64;
  eps_mmpde    = 1.e-3;
  smoothing    = smooth_opt(1);
  nsmooth      = 2;
  ksmooth      = 1.2e-1;
  avrg         = avrg_opt(1); 
  numtimesteps = 8;
  timestep     = 0.005;
  time 	       = time_start;
  eps_omega    = 1.e-2;
  interp_method= interp_opt(1);
  
  %Computational mesh is always uniform; ksi e [0,1]
  comp_mesh = linspace( 0, 1, nnodes)';
  
  %Initialize physical mesh as uniform
  physical_mesh = linspace( xa, xb, nnodes )';

  %Boundary and initial conditions
  ua = 0.0;
  ub = 0.0;
  u0 = exact_sol(physical_mesh,time,velocity);
  u1 = u0;

  %PDE coefficients arrays
  vel = velocity*ones(nnodes,1);
  dif = diffusion*ones(nnodes,1);

  %Time stepping
  for ts = 1:numtimesteps
    init_u0 = u0;
    init_mesh = physical_mesh;
    for niter = 0:niter_coupl
      %Solve nonlinear MMPDE to get final mesh
      [physical_mesh,u0,u1,omega,l2res,itr_mmpde] = adapt_mesh(nnodes,xa,xb,comp_mesh,init_u0,u1,init_mesh,physical_mesh,niter_mmpde,eps_mmpde,eps_omega,avrg,smoothing,nsmooth,ksmooth,interp_method);
      %Solve PDE using adaptive mesh
      [u1] = solve_ADE(nnodes,vel,dif,physical_mesh,timestep,u0,ua,ub,"arithmetic");
    end
    %u0 = interp1(init_mesh,init_u0,physical_mesh,interp_method,'extrap');
    %[u1] = solve_ADE(nnodes,vel,dif,physical_mesh,timestep,u0,ua,ub,"arithmetic");
    u0 = u1;
    time = time+timestep;
    fprintf('time = %f | l2res_mmpde = %e | iter_mmpde = %i\n',time,l2res,itr_mmpde);
  end

  %CV dx
  x = physical_mesh;
  delx = zeros(nnodes,1);
  delx(1,1) = x(2)-x(1);
  for i = 2:nnodes-1
    delx(i,1) = 0.5*(x(i+1)-x(i-1));
  end
  delx(nnodes,1) = x(nnodes)-x(nnodes-1);
  delx = delx/max(delx);


  %"Continuous" variables (for plot)
  dense_mesh    = linspace(xa,xb,5000)';
  u_exact       = exact_sol(dense_mesh,time,velocity);
  u_initial     = exact_sol(dense_mesh,time_start,velocity);
  grad_dense    = gradient(u_exact,dense_mesh);
  [omega_dense] = monitor_fun(grad_dense,eps_omega); 
  %Plot
  figure()
  plot(dense_mesh, u_initial, '-')
  hold on
  plot(dense_mesh, u_exact, '-')
  hold on
  plot(physical_mesh, u1, '-s')
  hold on
  plot(physical_mesh,comp_mesh, ':')
  hold on
  plot(physical_mesh, delx, '--')
  hold on
  plot(physical_mesh, 0.5, 'd')
  %hold on
  %plot(physical_mesh, grad_u0, '-s')
  xlim([xa,xb]);
  ylim([-0.5,1.5]);
  legend('u_{initial}','{u}_{exact}','{u}_{num}','{\xi}(x_j)','{\Delta x}','x_j');
  title('PDE')

  %%Plot
  %if(niter_mmpde~=0)
  %  figure()
  %  plot(dense_mesh, omega_dense, '-')
  %  hold on
  %  plot(physical_mesh, omega, '-s')
  %  hold on
  %  plot(physical_mesh,comp_mesh, '-')
  %  hold on
  %  plot(physical_mesh, 0.5, 'd')
  %  xlim([xa,xb]);
  %  legend('{\omega}_{exact}','{\omega}(x_j)','{\xi}(x_j)','x_j');
  %  title('MMPDE')
  %end

