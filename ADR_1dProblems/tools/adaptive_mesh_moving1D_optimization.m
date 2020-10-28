function [profile, final_mesh] = adaptive_mesh_moving1D_optimization(physical_grid, diff, diff_x, advection, advection_x, reaction, reaction_x, rhs, rhs_x )

    a = @( x ) diff;
    aprime = @( x ) diff_x * x;
    b = @( x ) advection + advection_x * x; 
    c = @( x ) reaction + reaction_x * x;
    f = @( x ) rhs + rhs_x + x; 
    
    [A,rhs] = prob_assemble( a, aprime, b, c, f, physical_grid );
    profile = A\rhs;    

    gradient_solution = Centered_gradient(physical_grid, profile);

    omega = monitor_function1D(physical_grid, gradient_solution); 
    omega = [1.0; omega; 1.0];
%     for i = 1:10
%         omega = smoothdata(omega, 'gaussian', kernel_width);
%     end     
%          

    [A_adaptive,rhs_adaptive] = elliptic_mesh_moving_problem_FINITE_VOLUMES( omega, physical_grid );
    map = A_adaptive\rhs_adaptive;    

    final_mesh = map; 

return 
end
