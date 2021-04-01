%L2 norm of discretized PDE residuals 
function [l2res] = residuals(x,A,b)

  l2res = 0.0;
  n = size(x,1); 
  
  %Loop over internal CVs (Dirichlet B.C. should be satisfied exactly) 
  for i = 2 : n-1
    res = A(i,i+1)*x(i+1) + A(i,i)*x(i) + A(i,i-1)*x(i-1) - b(i);
    l2res = l2res + res*res;
  end

  l2res = sqrt(l2res)/n;
