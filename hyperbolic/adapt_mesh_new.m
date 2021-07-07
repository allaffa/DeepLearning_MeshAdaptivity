function [adapt_mesh,l2res,iter] = adapt_mesh_new(bc_method, nnodes,xa,xb,comp_mesh,u0,uCurr,init_mesh,mesh,niter,eps_nonlin,eps_omega,avrg,smoothing,nsmooth,ksmooth,interp_method,dtau,Mvk,iur)
  dtau_init = dtau;

  %Initialize
  adapt_mesh  = mesh;
  init_uCurr  = uCurr;
  omega = zeros(nnodes,1);

  %%%Smooth profile function
  %for i=1:4
  %  if smoothing=="weightavrg" uCurr(:,i) = smooth_fun(uCurr(:,i),nsmooth,1);, end
  %  if smoothing=="gaussian"   uCurr(:,i) = smoothdata(uCurr(:,i), 'gaussian', ceil(ksmooth*length(adapt_mesh)));,end
  %end

  %Solve nonlinear MMPDE to get final mesh
  for iter = 1:niter
    %Calculate omega using mesh from previous iteration
    omega   = monitor_fun_new(uCurr,eps_omega,adapt_mesh);
    %%Smooth monitor function
    if smoothing=="weightavrg" omega = smooth_fun(omega,nsmooth,1); end
    if smoothing=="gaussian"   omega = smoothdata(omega, 'gaussian', ceil(ksmooth*length(adapt_mesh)));end
    
    % boundary condition
    if( isequal( bc_method, 'periodic' ) )
      omega(length(omega)) = omega(1);
    end
    
    %Solve MMPDE
    mesh = adapt_mesh;
    [adapt_mesh,A,b] = solve_MMPDE_new(omega,comp_mesh,nnodes,xa,xb,avrg,dtau,mesh,iur);
    
    if( Mvk == 1)
      %%check substep size time to prevent mesh tangling
      %[dtau_var] = substep_dtnk(adapt_mesh, mesh, dtau)
      %while( dtau > dtau_var )
      %  dtau = dtau/2.0;

      %  % deny adapt_mesh and continues
      %  [adapt_mesh,A,b] = solve_MMPDE(omega,comp_mesh,nnodes,xa,xb,avrg,dtau,mesh,iur);
        
      %  %check substep size time to prevent mesh tangling
      %  [dtau_var] = substep_dtnk(adapt_mesh, mesh, dtau);
      %end
      %dtau = dtau_init;
      
      %check substep size time to prevent mesh tangling
      [flag] = substep_dtnk_new(adapt_mesh, mesh);
      while( flag < 1 )
        dtau = dtau/2.0;

        % deny adapt_mesh and continues
        [adapt_mesh,A,b] = solve_MMPDE_new(omega,comp_mesh,nnodes,xa,xb,avrg,dtau,mesh,iur);
        
        %check substep size time to prevent mesh tangling
        [flag] = substep_dtnk_new(adapt_mesh, mesh);
      end
      dtau = dtau_init;
    end
    
    %Interpolate results from old to new mesh
      %uCurr(:,i) = interp1(init_mesh,init_uCurr(:,i),adapt_mesh,interp_method,'extrap');
      uCurr = interp1(mesh,uCurr,adapt_mesh,interp_method,'extrap');
      %u0(:,i) = interp1(mesh,u0(:,i),adapt_mesh,interp_method,'extrap');
    
    %Check convergence of nonlinear residuals (new coeff. with old mesh)
    [l2res] = residuals(mesh,A,b);
    %Normalize residual with initial residual
%     if iter == 1 l2res_1 = l2res;, end;
%     l2res = l2res/(l2res_1+1.e-12); %add small num for stability (if l2res_1 ~ 0)
    if l2res<eps_nonlin break, end;
  end

