function [final_mesh] = adaptive_mesh_moving1D_optimization( nnodes, physical_grid, solution, kernel_width)
   
    profile = solution; 
    
    gradient_solution = Centered_gradient(physical_grid, profile);

    omega = monitor_function1D(physical_grid, gradient_solution); 
    omega = [1.0; omega; 1.0];
    for i = 1:1
        omega = smoothdata(omega, 'gaussian', kernel_width);
    end     
        
    uniform_physical_grid = linspace(physical_grid(1), physical_grid(end), nnodes)';
    omega_uniform_physical_grid = interp1q(physical_grid, omega, uniform_physical_grid);
    omega_uniform_physical_grid = omega_uniform_physical_grid';    

    [A_adaptive,rhs_adaptive] = elliptic_mesh_moving_problem_FINITE_VOLUMES( omega_uniform_physical_grid, uniform_physical_grid );
    map = A_adaptive\rhs_adaptive;    

    final_mesh = map; 

return 
end