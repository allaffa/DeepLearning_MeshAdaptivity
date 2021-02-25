%Solution of advection diffusion equation (ADE):
%du/dt = -d/dx(vu) + d/dx (D du/dx) = 0 
%Currently Dirichlet b.c.
%2nd order FVM on physical (non-uniformn) grid
%FV boundaries are defined as midpoint between nodes
function [u] = solve_ADE(n,v,D,x,dt,u0,ua,ub,avrg)

  A = zeros(n,n);
  b = zeros(n,1);

  %Fill matrix A and rhs b
  for i = 2 : n-1
    %Vertex centered
    x_e = 0.5*(x(i+1)+x(i)); 
    x_w = 0.5*(x(i-1)+x(i)); 
    Vp  =(x_e-x_w);
    [v_w,v_e] = calc_flux(v(i-1),v(i),v(i+1),avrg,0,0);
    [D_w,D_e] = calc_flux(D(i-1),D(i),D(i+1),avrg,0,0);
    %Calculate discretized coefficients
    A(i,i+1) = + 0.5*v_e - D_e / (x(i+1)-x(i));
    A(i,i-1) = - 0.5*v_w - D_w / (x(i)-x(i-1));
    A(i,i)   = -( A(i,i+1) + A(i,i-1) );
    A(i,i)   = A(i,i) + Vp/dt; %temporal contibution
    %A(i,i)   = A(i,i) + 0.5*(-v_e+v_w); %compressible flow 
    b(i) = u0(i)*Vp/dt;
  end

  %Boundary conditions
  A(1,1) = 1.0; 
  A(n,n) = 1.0;
  b(1)   = ua;
  b(n)   = ub;

  %Solve linear system of equations
  u = A\b;
