function [final_mesh, uniform_grid, omega] = adaptive_mesh_moving1D_optimization( physical_grid, profile, kernel_width )
    
    %Map solution on uniform grid
    uniform_grid = linspace(physical_grid(1), physical_grid(end), length(physical_grid))';
    profile_uniform_grid = pchip(physical_grid, profile, uniform_grid);  
    for i = 1:1
        profile_uniform_grid = smoothdata( profile_uniform_grid, 'gaussian', ceil(kernel_width) );
    end         
    
%     gradient_profile = Centered_gradient(uniform_grid, profile_uniform_grid);
    gradient_profile = gradient(profile_uniform_grid, uniform_grid);
    
    omega = monitor_function1D(uniform_grid, gradient_profile); 
%     for i = 1:1
%         omega = smoothdata(omega, 'gaussian', ceil(kernel_width*length(uniform_grid)));
%     end     
%     omega = omega';
    
    [A_adaptive,rhs_adaptive] = elliptic_mesh_moving_problem_FINITE_VOLUMES( omega, uniform_grid);
    map = A_adaptive\rhs_adaptive; 
 

    final_mesh = map'; 


return 
end
