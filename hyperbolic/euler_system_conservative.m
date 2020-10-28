close all
clc

x0 = 0;
x1 = 1;
nnodes = 201;

t0 = 0.0;
tF = 0.1;
gamma = 1.4;

upwind = true;
verbose = true;
adaptive_mesh = true;

kernel_width = nnodes * 0.001;


physical_grid_uniform = linspace( x0, x1, nnodes )';

rho_initial = zeros(size(physical_grid_uniform));
u_initial = zeros(size(physical_grid_uniform));
p_initial = zeros(size(physical_grid_uniform));

for i = 1:length(physical_grid_uniform)
    
    if(physical_grid_uniform(i)<0.5) 
        rho_initial(i) = 8.0;
        u_initial(i) = 0.0;
        p_initial(i) = 10.0/gamma;
    else
        rho_initial(i) = 1.0;
        u_initial(i) = 0.0;
        p_initial(i) = 1/gamma;
    end
    
end

[physical_grid, u1_new, u2_new, u3_new, p_new, t_current, omega] = Euler_equations_solver_conservative(gamma, physical_grid_uniform, rho_initial, u_initial, p_initial, t0, tF, upwind, verbose, 1, ...
    kernel_width, adaptive_mesh);

figure()
plot(physical_grid, u1_new, '-o', 'Linewidth', 1)
hold on
plot(physical_grid, 4.0, '-o', 'Linewidth', 1)
xlim([physical_grid_uniform(1) physical_grid_uniform(end)])
title('Density')
set(gca, 'fontsize', 40)
title(strcat('Density - time = ', num2str(t_current)))

% figure()
% plot(physical_grid, u2_new./u1_new, '-o', 'Linewidth', 1)
% xlim([physical_grid_uniform(1) physical_grid_uniform(end)])
% title('Velocity')
% set(gca, 'fontsize', 40)
% title(strcat('Velocity - time = ', num2str(t_current)))
% 
% figure()
% plot(physical_grid, (u3_new - 1/2*u2_new.*u2_new./u1_new)./u1_new, '-o', 'Linewidth', 1)
% xlim([physical_grid_uniform(1) physical_grid_uniform(end)])
% title('Energy')
% set(gca, 'fontsize', 40)
% title(strcat('Energy - time = ', num2str(t_current)))
% 
% 
figure()
plot(physical_grid, p_new, '-o', 'Linewidth', 1)
xlim([physical_grid_uniform(1) physical_grid_uniform(end)])
title('Pressure')
set(gca, 'fontsize', 40)
title(strcat('Pressure - time = ', num2str(t_current)))


figure()
plot(physical_grid_uniform, physical_grid, '-o', 'Linewidth', 1)
xlim([physical_grid_uniform(1) physical_grid_uniform(end)])
ylim([physical_grid_uniform(1) physical_grid_uniform(end)])
title('New mesh')
set(gca, 'fontsize', 40)


figure()
plot(physical_grid_uniform, gradient_profile, '-o', 'Linewidth', 1)
xlim([physical_grid_uniform(1) physical_grid_uniform(end)])
title('Gradient_profile')
set(gca, 'fontsize', 40)


% figure()
% plot(physical_grid, 1.0, '-o', 'Linewidth', 1)
% xlim([physical_grid_uniform(1) physical_grid_uniform(end)])
% title('Flat new mesh plot')
% set(gca, 'fontsize', 40)

figure()
plot(physical_grid_uniform, omega, '-o', 'Linewidth', 1)
xlim([physical_grid_uniform(1) physical_grid_uniform(end)])
title('Monitor function')
set(gca, 'fontsize', 40)








