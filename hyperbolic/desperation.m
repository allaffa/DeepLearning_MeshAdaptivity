close all
clc

x0 = 0;
x1 = 1;
nnodes = 401;

t0 = 0.0;
n_time_steps = 10;
dt = 0.005;
t = t0;

kernel_width = nnodes * 1e-3;

% Analytical solution
sol_exact_fun = @(x, t, vel) 0.0 + 1.0 * ( ( (x >= 0.4 + vel * t) & (x <= 0.5 + vel * t ) ) );
% sol_exact_fun = @(x, t, vel) 0.0 + 1.0 * ( (x >= 0.1 + vel * t) );
% 
velocity = 1.0;

% time integration (Explicit Eulfer: alpha=0.0; Implicit Eulfer: alpha=1.0; Crank-Nicholson: alpha = 0.5)
alpha = 0.5;

physical_grid = linspace( x0, x1, nnodes )';

profile_exact = sol_exact_fun(physical_grid, t, velocity);
[final_mesh, uniform_grid, omega] = adaptive_mesh_moving1D_optimization( physical_grid, profile_exact, kernel_width );
figure()
set(gca, 'fontsize', 45);
hold on
plot(max(0,final_mesh), 0.5, 'o', 'Linewidth', 3)
% plot(max(0,physical_grid_uniform), omega, '-o', 'Linewidth', 1)
%ylim([-0.2 1.2])
xlabel('Coordinates in the physical domain')
ylabel('Solution values')
plot(max(0,physical_grid), profile_exact, '-', 'linewidth', 1);
title(strcat('Time: ', num2str(t)))
% plot(max(0,uniform_grid), omega)







