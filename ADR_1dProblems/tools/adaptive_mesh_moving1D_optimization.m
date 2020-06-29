function [initial_sol, original_gradient, final_mesh] = adaptive_mesh_moving1D_optimization(grid, epsilon, max_iter, diff, diff_x, advection, advection_x, reaction, reaction_x, rhs, rhs_x, verbose)

    a = @( x ) diff;
    aprime = @( x ) diff_x * x;
    b = @( x ) advection + advection_x * x; 
    c = @( x ) reaction + reaction_x * x;
    f = @( x ) rhs + rhs_x + x; 
    
    % n0 = initial value for the number of mesh nodes
    computational_domain = grid;
    x0 = grid(1);
    x1 = grid(end);
    
    original_grid = grid;
    discrepancy = 1;
    count_iter = 0;
    
    physical_domain = computational_domain; 
    
    while( discrepancy>epsilon && count_iter<max_iter )    
           
        [A,rhs] = prob_assemble( a, aprime, b, c, f, physical_domain );
        sol = A\rhs;    

        discrepancy = 1;
        initial_sol = sol;

        gradient_solution = Centered_gradient(physical_domain, sol);
        original_gradient = gradient_solution';

        omega = monitor_function1D(physical_domain, gradient_solution); 
        omega = [1.0; omega; 1.0];
        for i = 1:100
            omega = smoothdata(omega, 'gaussian', 5);
        end     
        omega = omega';
        gradient_omega = zeros(size(omega));
        gradient_omega = Centered_gradient(physical_domain, omega');
        gradient_omega = [0.0; gradient_omega; 0.0];
        gradient_omega = gradient_omega';
%         gradient_omega = 1e-4 * gradient_omega;
        [A_adaptive,rhs_adaptive] = elliptic_mesh_moving_problem( omega, gradient_omega, physical_domain' );
        map = A_adaptive\rhs_adaptive;
        
%         for i = 1:1
%             map = smoothdata(map, 'gaussian', 5);
%         end    
        
        discrepancy = max(abs(physical_domain - map));

        physical_domain = map; 
        count_iter = count_iter+1;
        
    end
                
    final_mesh = physical_domain'; 

return 
end
