  function [initial_sol, original_gradient, final_mesh] = adaptive_mesh_moving1D(grid, epsilon, diff, diff_x, advection, advection_x, reaction, reaction_x, rhs, rhs_x, verbose)

    count_iter = 0;
    
    a = @( x ) diff;
    aprime = @( x ) diff_x * x;
    b = @( x ) advection + advection_x * x; 
    c = @( x ) reaction + reaction_x * x;
    f = @( x ) rhs + rhs_x + x; 
    
    % n0 = initial value for the number of mesh nodes
    computational_domain = grid;
    x0 = grid(1);
    x1 = grid(end);
    
    physical_domain = computational_domain;    % Construct and Solve the problem on the coarse mesh
    [A,rhs] = prob_assemble( a, aprime, b, c, f, physical_domain );
    sol = A\rhs;
    original_gradient = Centered_gradient(physical_domain, sol);
    
    discrepancy = 1;
    initial_sol = sol;
    
    while( discrepancy>epsilon )
        
        gradient_solution = Centered_gradient(physical_domain, sol);
        omega = monitor_function1D(physical_domain, gradient_solution); 
        
        temp_physical_domain = zeros(size(computational_domain,1),1);
        temp_physical_domain(1) = x0;
        
        for i=1:size(omega,1)          
            temp_physical_domain(i+1) = temp_physical_domain(i) + 1/omega(i) ;     
        end
        
        temp_physical_domain(end) = temp_physical_domain(end-1) + 1/omega(end-1);
        
        % renormalize 
        temp_physical_domain = x0 + temp_physical_domain / (temp_physical_domain(end) - temp_physical_domain(1)) * (x1-x0);
        
        
        discrepancy = max( physical_domain - temp_physical_domain ); 
        physical_domain = temp_physical_domain; 
        
        [A,rhs] = prob_assemble( a, aprime, b, c, f, physical_domain );
        sol = A\rhs;
        
        count_iter = count_iter + 1; 
        
        if(verbose==true)
            display(strcat('Iteration counter: ', num2str(count_iter), ' - discrepancy: ', num2str(discrepancy)));
        end
         
    end
    
    final_mesh = physical_domain; 

return 
end