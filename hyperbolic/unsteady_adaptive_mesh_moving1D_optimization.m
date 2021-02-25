function [final_mesh, uniform_grid, omega] = unsteady_adaptive_mesh_moving1D_optimization( physical_grid, profile, kernel_width, dt )
    
    %Map solution on uniform grid
    uniform_grid = linspace(physical_grid(1), physical_grid(end), length(physical_grid))';
    profile_uniform_grid = interp1(physical_grid, profile, uniform_grid);  
    gradient_profile = gradient(profile_uniform_grid, uniform_grid);
    
    omega = monitor_function1D(uniform_grid, gradient_profile); 
    for i = 1:1
        omega = smoothdata(omega, 'gaussian', ceil(kernel_width));
    end     
    omega = omega';
    
    dxsi = min(diff(uniform_grid));
    

    final_mesh = zeros( size(physical_grid) ); 
    final_mesh(1) = physical_grid(1);
    final_mesh(end) = physical_grid(end);
    
    for time_step=1:10
 
        for i = 2:length(physical_grid)-1

            omega_e = ( omega(i+1)+omega(i) )/2;
            omega_w = ( omega(i)+omega(i-1) )/2;

            final_mesh(i) = physical_grid(i) + dt * ( omega_e * (physical_grid(i+1)-physical_grid(i)) -  omega_w * (physical_grid(i)-physical_grid(i-1)) )/(dxsi*dxsi); 
        end
        
        physical_grid = final_mesh;
    
    end


return 
end
