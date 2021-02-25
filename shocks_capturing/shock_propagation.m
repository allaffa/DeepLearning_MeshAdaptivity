pyExec	= '/opt/anaconda3/bin/python3';
[ver, exec, loaded]	= pyversion(pyExec);

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

%Problem definiton
  %exact_sol = @(x, t, v) 0.0 + 1.0 * ( ( (x >= 0.25 + v * t) & (x <= 0.45 + v * t ) ) ) + 1.0 * ( ( (x >= 0.65 + v * t) & (x <= 0.85 + v * t ) ) );
  exact_sol = @(x, t, v) 0.0 + 1.0 * ( ( (x >= 0.30 + v * t) & (x <= 0.50 + v * t ) ) ) ;
  
  %Physical domain boundaries
  xa = 0.0;
  xb = 1.0;

  %Parameters
  smooth_opt   = ["none","weightavrg","gaussian"];
  avrg_opt     = ["arithmetic","harmonic","geometric"];
  interp_opt   = ["linear","pchip"];
  nnodes       = 201;
  time_start   = 0.0;
  velocity     = 1.e-2;
  diffusion    = 0.e-5;
  niter_mmpde  = 4;
  niter_coupl  = 2;
  smoothing    = smooth_opt(1);
  nsmooth      = 2;
  ksmooth      = 1.2e-1;
  avrg         = avrg_opt(1); 
  numtimesteps = 16;
  timestep     = 0.5;
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
  
  AI_method = 'DL';

  %Time stepping
  for ts = 1:numtimesteps
      
    display(strcat("Time step: ",num2str(ts)," of ",num2str(numtimesteps)))
    
    init_mesh = physical_mesh;
      
    if strcmp(AI_method, 'DL')
        physical_mesh = double(py.array.array('d',py.numpy.nditer(python_module.evaluate_model_load(u0))));
        physical_mesh = sort(physical_mesh);
    else
        standard_deviation =  1e-2;
        physical_mesh = importance_sampling_mesh(physical_mesh, u0, standard_deviation);
    end

    u0 = interp1(init_mesh,u0,physical_mesh,'linear','extrap');

    %Solve PDE using adaptive mesh
    [u1] = solve_ADE(nnodes,vel,dif,physical_mesh,timestep,u0,ua,ub,"arithmetic");
    
    u0 = u1;
    time = time+timestep;
    
  end
  

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
  plot(physical_mesh, 0.5, 'd')
  %hold on
  %plot(physical_mesh, grad_u0, '-s')
  xlim([xa,xb]);
  ylim([-0.5,1.5]);
  legend('u_{initial}','{u}_{exact}','{u}_{num}','{\xi}(x_j)','x_j');
  title('PDE')

 