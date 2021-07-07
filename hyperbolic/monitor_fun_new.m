function [omega] = monitor_fun_new(f,eps,adapt_mesh)

  %omega = sqrt(1.0+eps*grad_f.*grad_f); 
%     grad_f(:,i) = gradient(f(:,i),adapt_mesh);
  grad_f = gradient(f,adapt_mesh);
  max_f  = max(abs(grad_f));
  
  n = length(grad_f);
  omega = zeros(n,1);

  for i = 1 : n
   omega(i,1) = sqrt(1.0+eps*grad_f(i,1)*grad_f(i,1));
%     omega(i,1) = sqrt( 1.0 + eps(1)*(grad_f(i,1)/(max_f(1)+1.e-12))^2 + eps(2)*(grad_f(i,2)/max_f(2))^2 + eps(3)*(grad_f(i,3)/max_f(3))^2 + eps(4)*(grad_f(i,4)/(max_f(4)+1.0e-12))^2);
  end

