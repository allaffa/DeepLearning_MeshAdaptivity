%Solution of linearized MMPDE: d/dksi (omega dx/dksi) = 0 
%2nd order FVM on (computational) uniform grid
function [final_mesh,A,b] = solve_MMPDE_new(omega,ksi,n,xa,xb,avrg,dtau,x0,iur)

  A = zeros(n,n);
  b = zeros(n,1);

  %Fill matrix A and rhs b
  for i = 2 : n-1
    %Average omega values
    [omega_w,omega_e] = calc_flux(omega(i-1),omega(i),omega(i+1),avrg,0,0);
    %Calculate discretized coefficients
    A(i,i+1) = omega_e / (ksi(i+1)-ksi(i));
    A(i,i-1) = omega_w / (ksi(i)-ksi(i-1));
    A(i,i)   = -( A(i,i+1) + A(i,i-1) +(ksi(i+1) - ksi(i-1))/(2*dtau) );
    %Rhs
    b(i) = -x0(i) * (ksi(i+1) - ksi(i-1))/(2*dtau);
  end

  %Under-relaxation
  for i = 2 : n-1
    b(i) = b(i) + (1.0-iur)/iur*A(i,i)*x0(i);
    A(i,i) = A(i,i)/iur;
  end

  %Boundary conditions
  A(1,1) = 1.0; 
  A(n,n) = 1.0;
  b(1)   = xa;
  b(n)   = xb;

  %Solve linear system of equations
  final_mesh = A\b;
