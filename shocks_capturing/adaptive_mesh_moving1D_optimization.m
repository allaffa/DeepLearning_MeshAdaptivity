function [adapted_physical_mesh] = adaptive_mesh_moving1D_optimization(physical_grid, profile, kernel_width)
    
    reference_grid = linspace(physical_grid(1), physical_grid(end), length(physical_grid));
    reference_grid = reference_grid';

    profile_reference_grid = profile;
        
    gradient_profile_reference_grid = Centered_gradient(reference_grid, profile_reference_grid);
    omega = monitor_function1D(reference_grid, gradient_profile_reference_grid); 
    omega = [1.0; omega; 1.0];
    
%     for i = 1:1
%         omega = smoothdata(omega, 'gaussian', kernel_width);
%     end     
%     omega = omega';
    
    [A_adaptive,rhs_adaptive] = elliptic_mesh_moving_problem_FINITE_VOLUMES( omega, reference_grid, physical_grid );
    map = A_adaptive\rhs_adaptive;

    adapted_physical_mesh = map; 

return 
end
