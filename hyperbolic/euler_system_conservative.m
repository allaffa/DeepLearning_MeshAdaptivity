close all
clear

%%

pyExec	= '/opt/anaconda3/bin/python3';
% [ver, exec, loaded]	= pyversion(pyExec);

python_path = py.sys.path;

% Add folders to python system path to loead version 3.7 of Python.
if count(python_path, pyExec) == 0
    insert(py.sys.path, int64(0), pyExec);
end

% Verify that the Python 3.7 version from Anaconda is successfully loaded
[ver2, exec2, loaded2]	= pyversion;
assert(loaded2==1);

python_module = py.importlib.import_module('model_load_script');

%%

x0 = 0;
x1 = 0.6;
nghost = 3;
nodes = 201;
nnodes = nodes + 2*nghost;
xc = 0.5;

t0 = 0.0;
tF = 1.0;
gamma = 1.4;

upwind = true;
verbose = true;
adaptive_mesh = "none";

kernel_width = nnodes * 0.05;

physical_grid_uniform = linspace( x0-nghost*(x1-x0)/(nodes-1), x1+nghost*(x1-x0)/(nodes-1), nnodes )';

rho_initial = zeros(size(physical_grid_uniform));
u_initial = zeros(size(physical_grid_uniform));
p_initial = zeros(size(physical_grid_uniform));

for i = 1:length(physical_grid_uniform)
    
%     %% sod
%     if(physical_grid_uniform(i)<0.5) 
%         rho_initial(i) = 1.0;
%         u_initial(i) = 0.0;
%         p_initial(i) = 1.0;
%     else
%         rho_initial(i) = 0.125;
%         u_initial(i) = 0.0;
%         p_initial(i) = 0.1;
%     end
    
%     %%  Woodward–Colella blast wave
%     if(physical_grid_uniform(i)<=0.1) 
%         rho_initial(i) = 1.0;
%         u_initial(i) = 0.0;
%         p_initial(i) = 1000.0;
%     elseif(physical_grid_uniform(i)>0.1 && physical_grid_uniform(i)<=0.9) 
%         rho_initial(i) = 1.0;
%         u_initial(i) = 0.0;
%         p_initial(i) = 0.01;
%     else
%         rho_initial(i) = 1.0;
%         u_initial(i) = 0.0;
%         p_initial(i) =100.0;
%     end

    %% sedov
    if( i >= 1 + nghost && i <= 2 + nghost )
%     if(physical_grid_uniform(i)<=1.0/(nodes-1) + eps && physical_grid_uniform(i) >=0 - eps ) 
        rho_initial(i) = 1.0;
        u_initial(i) = 0.0;
        cell_mass = rho_initial(i) * 0.6/(nodes-1) ;
        spe = 0.0673185/(2.0) / cell_mass;
        p_initial(i) = (gamma-1.0) * rho_initial(i) * spe;
    else
        rho_initial(i) = 1.0;
        u_initial(i) = 0.0;
%         cell_mass = rho_initial(i) * 1.0/(nodes-1) ;
%         spe = 0.0673185e-6 / cell_mass;
        p_initial(i) = 0.0;
    end
    
end

tic
[physical_grid, u1_new, u2_new, u3_new, p_new, t_current] = Euler_equations_solver_conservative_WENO5(gamma, physical_grid_uniform, rho_initial, u_initial, p_initial, t0, tF, nghost);
% [physical_grid, u1_new, u2_new, u3_new, p_new, t_current] = Euler_equations_solver_conservative(gamma, physical_grid_uniform, rho_initial, u_initial, p_initial, t0, tF, verbose, kernel_width, adaptive_mesh, python_module);
toc

%
adaptive_mesh = "none";
tic
[physical_grid2, u1_new2, u2_new2, u3_new2, p_new2, t_current2, time_history2, grid_history2] = Euler_equations_solver_conservative_WENO5_nonuniform(gamma, physical_grid_uniform, rho_initial, u_initial, p_initial, t0, tF, nghost, verbose, kernel_width, adaptive_mesh, python_module);
% [physical_grid2, u1_new2, u2_new2, u3_new2, p_new2, t_current2] = Euler_equations_solver_conservative(gamma, physical_grid_uniform, rho_initial, u_initial, p_initial, t0, tF, verbose, kernel_width, adaptive_mesh, python_module);
toc

%
adaptive_mesh = "DL";
tic
[physical_grid3, u1_new3, u2_new3, u3_new3, p_new3, t_current3, time_history3, grid_history3] = Euler_equations_solver_conservative_WENO5_nonuniform(gamma, physical_grid_uniform, rho_initial, u_initial, p_initial, t0, tF, nghost, verbose, kernel_width, adaptive_mesh, python_module);
% [physical_grid3, u1_new3, u2_new3, u3_new3, p_new3, t_current3] = Euler_equations_solver_conservative(gamma, physical_grid_uniform, rho_initial, u_initial, p_initial, t0, tF, verbose, kernel_width, adaptive_mesh, python_module);
toc

%%

%% Plots

% [data] = analytic_sod(t_current, nnodes);
[data] = exact_sedov_1d(gamma, linspace(0, 0.6, 241)');

figure()
plot(physical_grid_uniform, u1_new , 'b-o', 'Linewidth', 6)
hold on
plot(physical_grid3, u1_new3, 'g-o', 'Linewidth', 6)
plot(physical_grid2, u1_new2, 'r-o', 'Linewidth', 6)
plot(data.c_x, data.rho, 'k', 'Linewidth', 6)
xlim([physical_grid_uniform(1) physical_grid_uniform(end)])
title('Density')
set(gca, 'fontsize', 72)
h = legend('Uniform mesh', 'DL adaptivity', 'Standard adaptivity', 'Exact solution');
set(h, 'Fontsize', 50)
% title(strcat('Density - time = ', num2str(t_current)))
title(strcat('Density - time = ', num2str(1.0)))


figure()
plot(physical_grid_uniform, u2_new./u1_new, 'b-o', 'Linewidth', 6)
hold on
plot(physical_grid3, u2_new3./u1_new3, 'g-o', 'Linewidth', 6)
plot(physical_grid2, u2_new2./u1_new2, 'r-o', 'Linewidth', 6)
plot(data.x, data.u, 'k', 'Linewidth', 6)
xlim([physical_grid_uniform(1) physical_grid_uniform(end)])
title('Velocity')
set(gca, 'fontsize', 72)
h = legend('Uniform mesh', 'DL adaptivity', 'Standard adaptivity', 'Exact solution');
set(h, 'Fontsize', 50)
% title(strcat('Velocity - time = ', num2str(t_current)))
title(strcat('Velocity - time = ', num2str(1.0)))

figure()
plot(physical_grid_uniform, (u3_new - 1/2*u2_new.*u2_new./u1_new)./u1_new, 'b-o', 'Linewidth', 6)
hold on
plot(physical_grid3, (u3_new3 - 1/2*u2_new3.*u2_new3./u1_new3)./u1_new3, 'g-o', 'Linewidth', 6)
plot(physical_grid2, (u3_new2 - 1/2*u2_new2.*u2_new2./u1_new2)./u1_new2, 'r-o', 'Linewidth', 6)
plot(data.c_x, data.e, 'k', 'Linewidth', 6)
xlim([physical_grid_uniform(1) physical_grid_uniform(end)])
title('Energy')
set(gca, 'fontsize', 72)
h = legend('Uniform mesh', 'DL adaptivity', 'Standard adaptivity', 'Exact solution');
set(h, 'Fontsize', 50)
% title(strcat('Energy - time = ', num2str(t_current)))
title(strcat('Energy - time = ', num2str(1.0)))

% 
figure()
plot(physical_grid_uniform, p_new, 'b-o', 'Linewidth', 6)
hold on
plot(physical_grid3, p_new3, 'g-o', 'Linewidth', 6)
plot(physical_grid2, p_new2, 'r-o', 'Linewidth', 6)
plot(data.c_x, data.P, 'k', 'Linewidth', 6)
xlim([physical_grid_uniform(1) physical_grid_uniform(end)])
title('Pressure')
set(gca, 'fontsize', 72)
h = legend('Uniform mesh', 'DL adaptivity', 'Standard adaptivity', 'Exact solution');
set(h, 'Fontsize', 50)
% title(strcat('Pressure - time = ', num2str(t_current)))
title(strcat('Pressure - time = ', num2str(1.0)))


figure()
plot(physical_grid_uniform, physical_grid2, '-o', 'Linewidth', 1)
xlim([physical_grid_uniform(1) physical_grid_uniform(end)])
ylim([physical_grid_uniform(1) physical_grid_uniform(end)])
title('New mesh')
set(gca, 'fontsize', 72)

figure
plot(diff(physical_grid2), '-o', 'Linewidth', 6)
hold on
plot(diff(physical_grid3), '-o', 'Linewidth', 6)
legend('Standard adaptivity', 'DL adaptivity')
title('Mesh spacing')
set(gca, 'fontsize', 72)

figure()
for i=1:size(grid_history2, 1)
    if mod(i,25)==0
        plot(grid_history2(i,:), time_history2(i), 'b.', 'Linewidth', 1);
        xlim([x0, x1])
        set(gca, 'Fontsize', 48)
        xlabel('Non-uniform mesh coordinates')
        ylabel('time step')
        title('Standard mesh adaptivity')
        hold on
    end
end

figure()
for i=1:size(grid_history3, 1)
    if mod(i,10)==0
        plot(grid_history3(i,:), time_history3(i), 'b.', 'Linewidth', 1);
        xlim([x0, x1])
        set(gca, 'Fontsize', 48)
        xlabel('Non-uniform mesh coordinates')
        ylabel('time step')
        title('DL adaptivity')
        hold on
    end
end





