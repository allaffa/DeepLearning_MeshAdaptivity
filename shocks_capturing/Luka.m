
n=10;
xsi = linspace(0,1,n); 
omega = ones(size(xsi));

gaussian =@(x) exp(-((x-0.2)/0.2).^2);

for i = 1:length(xsi)
    omega(i) = gaussian(xsi(i));
end

x = zeros(size(xsi)); 
x(n) = 5;
x(1) = 2;

A = zeros ( n, n );
rhs = zeros ( n, 1 );
%
%  The first equation is the left boundary condition, U(X(1)) = 0.0;
%
  A(1,1) = 1.0;
  rhs(1) = x(1);
%
%  Now gather the multipliers of U(I-1) to get the matrix entry A(I,I-1),
%  and so on.
%
  for i = 2 : n - 1

    A(i,i-1) = - ( omega(i)+omega(i-1) )/( 2*(xsi(i)-xsi(i-1)) );

    A(i,i+1) = - ( omega(i)+omega(i+1) )/( 2*(xsi(i+1)-xsi(i)) );
    
    A(i,i) = - ( A(i,i-1) + A(i,i+1) );

  end
%
%  The last equation is the right boundary condition, U(X(N)) = 0.0;
%
  A(n,n) = 1.0;
  rhs(n) = x(n);
%
%  Solve the linear system.
%
 
x = A\rhs;

new_gaussian =@(x) exp(-((x-2.2)/0.1).^2);
new_omega = ones(size(x));

for i = 1:length(xsi)
    new_omega(i) = new_gaussian(x(i));
end

x_uniform = linspace(x(1), x(end), length(x));
x_uniform = x_uniform';

new_omega_interp = interp1q(x, new_omega, x_uniform);

[A,rhs] = elliptic_mesh_moving_problem_FINITE_VOLUMES ( new_omega_interp, xsi, x_uniform );

x_new = A\rhs;


