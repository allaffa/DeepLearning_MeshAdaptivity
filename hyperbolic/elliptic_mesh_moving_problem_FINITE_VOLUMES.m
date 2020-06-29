function [A,rhs] = elliptic_mesh_moving_problem_FINITE_VOLUMES ( a, x )

%*****************************************************************************80
%
%% FD1D_BVP solves a two point boundary value problem.
%
%  Discussion:
%
%    The program uses the finite volume method to solve a BVP
%    (boundary value problem) in one dimension.
%
%    The problem is defined on the region X(1) <= x <= X(N).
%
%    The following differential equation is imposed in the region:
%
%      - d/dx a(x) du/dx = f(x)
%
%    where a(x), c(x), and f(x) are given functions.  We write out
%    the equation in full as
%
%      - a(x) * u''(x) - a'(x) * u'(x) + c(x) * u(x) = f(x)
%
%    At the boundaries, the following conditions are applied:
%
%      u(X(1)) = u_1
%      u(X(N)) = u_N
%
%    A set of N nodes is defined on this
%    interval, with X(1) < X(2) < ... < X(N).
%
%    We replace the function U(X) by a vector of values U(1)
%    through U(N), associated with the nodes.
%
%    The values of U(1) and U(N) are determined by the boundary conditions.
%
%    At each interior node I, we write an equation to help us determine
%    U(I).  We do this by approximating the derivatives of U(X) by
%    finite volumes.  Let us write XL, XM, and XR for X(I-1), X(I) and X(I+1).
%    Similarly we have UL, UM, and UR.  Other quantities to be evaluated at
%    X(I) = XM will also be labeled with an M:
%
%    - a_L * u_L + a_M * u_M - a_R * u_R = f_M 
%
%    These N-2 linear equations for the unknown coefficients complete the
%    linear system and allow us to compute the finite volume approximation
%    to the solution of the BVP.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    25 August 2010
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer N, the number of nodes.
%
%    Input, function A ( X ), evaluates a(x);
%
%    Input, function F ( X ), evaluates f(x);
%
%    Input, real X(N), the mesh points, which may be nonuniformly spaced.
%
%    Output, real U(N), the value of the finite difference approximation
%    to the solution.
%


%  generate a uniform grid fro the reference domain
  n = size(x,1);
  xsi = linspace(0, 1, n)';
  
% N.B. : the gradients are computed on the physical domain, NOT on the
% reference domain

%
%  Make room for the matrix A and right hand side b.
%
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

    A(i,i-1) = - ( a(i)+a(i-1) )/( 2*(xsi(i)-xsi(i-1)) );

    A(i,i+1) = - ( a(i)+a(i+1) )/( 2*(xsi(i+1)-xsi(i)) );
    
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
  return
end
