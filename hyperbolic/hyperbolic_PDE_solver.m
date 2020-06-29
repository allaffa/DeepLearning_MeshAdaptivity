function [sol] = hyperbolic_PDE_solver(x, sol, dt, c, alpha)

      n = length(x);
      A = zeros ( n, n );
      rhs = zeros ( n, 1 );
    %
    %  The first equation is the left boundary condition, U(X(1)) = 0.0;
    %
      A(1,1) = 1.0;
    %
    %  Now gather the multipliers of U(I-1) to get the matrix entry A(I,I-1),
    %  and so on.
    %
      for i = 2 : n - 1
          
        x_e = ( x(i)+x(i+1) )/2.0;
        
        x_w = ( x(i-1)+x(i) )/2.0;
        
        g_e = ( x_e -  x(i) )/( x(i+1) - x(i) );
        g_w = ( x(i) -  x_w )/( x(i) - x(i-1) );

        A(i,i-1) = - alpha * g_w;

        A(i,i+1) = + alpha * g_e;

        A(i,i) = ( x_e - x_w )/( c * dt ) + alpha * [ ( 1 - g_e ) - ( 1-g_w ) ];

        rhs(i) = [ (x_e - x_w)/(c * dt) - (1-alpha)*(1-g_e) + (1-alpha)*(1-g_w) ] * sol(i) + ...
                 [ - (1-alpha) * g_e ] * sol(i+1) + [ (1-alpha) * g_w ] * sol(i-1);
        
      end
    %
    %  The last equation is the right boundary condition, U(X(N)) = 0.0;
    %
      A(n,n) = 1.0;
    %
    %  Solve the linear system.
    %
    
    sol = A\rhs;

return