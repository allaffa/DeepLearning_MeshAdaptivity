close all
clear

x0 = 0;
x1 = 0.6;
nghost = 3;
nodes = 2401;
nnodes = nodes + 2*nghost;
xc = 0.5;

t0 = 0.0;
tF = 1.0;
gamma = 1.4;

upwind = true;
verbose = true;
adaptive_mesh = false;

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


[physical_grid, u1_new, u2_new, u3_new, p_new, t_current] = Euler_equations_solver_conservative_WENO5(gamma, physical_grid_uniform, rho_initial, u_initial, p_initial, t0, tF, nghost);
% [physical_grid, u1_new, u2_new, u3_new, p_new, t_current] = Euler_equations_solver_conservative(gamma, physical_grid_uniform, rho_initial, u_initial, p_initial, t0, tF, verbose, kernel_width, adaptive_mesh);

                                                                                               
[data] = analytic_sod(t_current, nnodes);

figure()
% plot(physical_grid, data.rho, 'Linewidth', 1)
hold on
%plot(exactsol.x, exactsol.rho)
plot(physical_grid, u1_new, 'r.', 'Linewidth', 2)
xlim([physical_grid_uniform(1) physical_grid_uniform(end)])
title('Density')
% set(gca, 'fontsize', 40)
% title(strcat('Density - time = ', num2str(t_current)))

figure()
% plot(physical_grid, data.u, 'Linewidth', 1)
hold on
plot(physical_grid, u2_new./u1_new, 'r.', 'Linewidth', 2)
xlim([physical_grid_uniform(1) physical_grid_uniform(end)])
title('Velocity')
% set(gca, 'fontsize', 40)
% title(strcat('Velocity - time = ', num2str(t_current)))

figure()
% plot(physical_grid, data.e, 'Linewidth', 1)
hold on
plot(physical_grid, (u3_new - 1/2*u2_new.*u2_new./u1_new)./u1_new, 'r.', 'Linewidth', 2)
xlim([physical_grid_uniform(1) physical_grid_uniform(end)])
title('Energy')
% set(gca, 'fontsize', 40)
% title(strcat('Energy - time = ', num2str(t_current)))

% 
figure()
% plot(physical_grid, data.P, 'Linewidth', 1)
hold on
plot(physical_grid, p_new, 'r.', 'Linewidth', 2)
xlim([physical_grid_uniform(1) physical_grid_uniform(end)])
title('Pressure')
% set(gca, 'fontsize', 40)
% title(strcat('Pressure - time = ', num2str(t_current)))


figure()
plot(physical_grid_uniform, physical_grid, '-o', 'Linewidth', 1)
xlim([physical_grid_uniform(1) physical_grid_uniform(end)])
ylim([physical_grid_uniform(1) physical_grid_uniform(end)])
title('New mesh')
set(gca, 'fontsize', 40)


% figure()
% plot(physical_grid, 1.0, '-o', 'Linewidth', 1)
% xlim([physical_grid_uniform(1) physical_grid_uniform(end)])
% title('Flat new mesh plot')
% set(gca, 'fontsize', 40)

% figure()
% plot(physical_grid_uniform, omega, '-o', 'Linewidth', 1)
% xlim([physical_grid_uniform(1) physical_grid_uniform(end)])
% title('Monitor function')
% set(gca, 'fontsize', 40)

save('sedov_weno5_2Kcell');
% save('sedov_weno5');





