function [omega] = monitor_fun(grad_f,eps)

  %omega = sqrt(1.0+eps*grad_f.*grad_f); 
  
  n = size(grad_f,1);
  omega = zeros(n,1);

  for i = 1 : n
    omega(i,1) = sqrt(1.0+eps*grad_f(i,1)*grad_f(i,1));
  end

