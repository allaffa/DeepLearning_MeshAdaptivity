function [final_mesh, uniform_grid, omega] = adaptive_mesh_moving1D_optimization( physical_grid, profile, kernel_width )
    
    %Map solution on uniform grid
    uniform_grid = linspace(physical_grid(1), physical_grid(end), length(physical_grid))';
    profile_uniform_grid = interp1(physical_grid, profile, uniform_grid);  
    gradient_profile = gradient(profile_uniform_grid, uniform_grid);
    
    omega = monitor_function1D(uniform_grid, gradient_profile); 
    for i = 1:1
        omega = smoothdata(omega, 'gaussian', ceil(kernel_width));
    end     
    omega = omega';
    
    [A_adaptive,rhs_adaptive] = elliptic_mesh_moving_problem_FINITE_VOLUMES( omega, uniform_grid);
    map = A_adaptive\rhs_adaptive; 
 

    final_mesh = map'; 


return 
end
