close all
clc

x0 = 0;
x1 = 1;
nnodes = 201;

t0 = 0.0;
n_time_steps = 2;
dt = 1e-7;
t = t0;

kernel_width = nnodes * 1e-3;

% Analytical solution
a = 0.5;
f_a = 1;
b = 0.55;
f_b = 4;
c = 0.6;
f_c = 1;
d = 0.65;
f_d = 2;
sol_exact_fun = @(x, t, vel) 0.0 + ( ( (x >= a + vel * t) & (x <= b + vel * t ) ) ) .* ( f_a + (f_b-f_a)/(b-a)*(x-a) ) ...
                + ( ( (x >= c + vel * t) & (x <= d + vel * t ) ) ) .* ( f_c + (f_d-f_c)/(d-c)*(x-c) );
               
% sol_exact_fun = @(x, t, vel) 0.0 + 1.0 * ( ( (x >= a + vel * t) & (x <= b + vel * t ) ) ) + 1.0 * ( ( (x >= c + vel * t) & (x <= d + vel * t ) ) );
% sol_exact_fun = @(x, t, vel) 0.0 + 1 * ( ( (x >= 0.01+ vel * t) ) );
% 
velocity = 1.0;

% time integration (Explicit Eulfer: alpha=0.0; Implicit Eulfer: alpha=1.0; Crank-Nicholson: alpha = 0.5)
alpha = 0.5;

physical_grid = linspace( x0, x1, nnodes )';

profile_exact = sol_exact_fun(physical_grid, t, velocity);
% [final_mesh, uniform_grid, omega] = adaptive_mesh_moving1D_optimization( physical_grid, profile_exact, kernel_width );
[final_mesh, uniform_grid, omega] = unsteady_adaptive_mesh_moving1D_optimization( physical_grid, profile_exact, kernel_width, dt );

figure()
set(gca, 'fontsize', 45);
hold on
plot(max(0,final_mesh), 0.5, 'o', 'Linewidth', 3)
% plot(max(0,physical_grid_uniform), omega, '-o', 'Linewidth', 1)
%ylim([-0.2 1.2])
xlabel('Coordinates in the physical domain')
ylabel('Solution values')
plot(max(0,physical_grid), profile_exact, '-', 'linewidth', 5);
title(strcat('Time: ', num2str(t)))
% plot(max(0,uniform_grid), omega)







