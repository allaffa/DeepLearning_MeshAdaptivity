function [adapt_mesh,u0,uCurr,omega] = adapt_mesh(nnodes,xa,xb,comp_mesh,u0,uCurr,mesh,niter,eps_omega,avrg,smoothing,nsmooth,ksmooth,interp_method)

  %Initialize
  adapt_mesh  = mesh;
  init_mesh   = mesh;
  init_uCurr  = uCurr;
  omega = zeros(nnodes,1);

  %Solve nonlinear MMPDE to get final mesh
  for iter = 1:niter
    %Calculate omega using mesh from previous iteration
    grad_uCurr = gradient(uCurr,adapt_mesh);
    omega   = monitor_fun(grad_uCurr,eps_omega);
    %%Smooth monitor function
    if smoothing=="weightavrg" omega = smooth_fun(omega,nsmooth,1); end
    if smoothing=="gaussian"   omega = smoothdata(omega, 'gaussian', ceil(ksmooth*length(adapt_mesh))); end
    %Solve MMPDE
    mesh = adapt_mesh;
    [adapt_mesh] = solve_MMPDE(omega,comp_mesh,nnodes,xa,xb,avrg);
    %Interpolate results from old to new mesh
    uCurr = interp1(init_mesh,init_uCurr,adapt_mesh,interp_method,'extrap');
    %uCurr = interp1(mesh,uCurr,adapt_mesh,interp_method,'extrap');
    %u0 = interp1(mesh,u0,adapt_mesh,interp_method,'extrap');
  end

  %Interpolate u0 solution on final mesh
  u0 = interp1(init_mesh,u0,adapt_mesh,interp_method,'extrap');

