%% Offline creation of the dataset

addpath('./tools')

num_experiments = 50000;

n0 = 201;
x0 = 0; 
x1 = 1;
epsilon = 1e-8;

physical_grid = linspace(x0, x1, n0)'; 

diffusion_coeff_range = [0 1];
advection_coeff_range = [-10 10];
reaction_coeff_range = [-1 1];
rhs_range = [-10 10];


% gradients = zeros(num_experiments, n0-1);
% initial_meshes = zeros(num_experiments, n0);
adapted_meshes = zeros(num_experiments, n0);
solutions = zeros(num_experiments, n0);
% initial_deltas = zeros(num_experiments, n0-1);
% final_deltas = zeros(num_experiments, n0-1);


verbose = true;

reference_grid = linspace(0,1,n0)';

for iter = 1:num_experiments
    
    diffusion = rand * (diffusion_coeff_range(2)-diffusion_coeff_range(1)) + diffusion_coeff_range(1); 
    diffusion_x = rand * (diffusion_coeff_range(2)-diffusion_coeff_range(1) + diffusion_coeff_range(1)); 
    advection = rand * (advection_coeff_range(2)-advection_coeff_range(1)) + advection_coeff_range(1); 
    advection_x = rand * (advection_coeff_range(2)-advection_coeff_range(1) + advection_coeff_range(1)); 
    reaction = rand * (reaction_coeff_range(2)-reaction_coeff_range(1)) + reaction_coeff_range(1); 
    reaction_x = rand * (reaction_coeff_range(2)-reaction_coeff_range(1) + reaction_coeff_range(1)); 
    rhs = rand * (rhs_range(2)-rhs_range(1)) + rhs_range(1);
    rhs_x = rand * (rhs_range(2)-rhs_range(1)) + rhs_range(1);
    
    display(strcat('Problem: ', num2str(iter), ' - of: ', num2str(num_experiments)));
    
%     [sol, original_gradient, final_mesh] = adaptive_mesh_moving1D(grid, epsilon, diff, diff_x, advection, advection_x, reaction, reaction_x, rhs, rhs_x, verbose);
    [sol, final_mesh] = adaptive_mesh_moving1D_optimization(physical_grid, diffusion, diffusion_x, advection, advection_x, reaction, reaction_x, rhs, rhs_x );
    
    adapted_meshes(iter,:) = final_mesh';
    solutions(iter,:) = sol'; 
%     initial_deltas(iter,:) = diff(grid)';
%     gradients(iter,:) = (diff(sol)')./initial_deltas(iter,:);
%     final_deltas(iter,:) = diff(final_mesh)';
    
end

% dlmwrite('data_set/initial_meshes.txt',initial_meshes,'delimiter','\t')
dlmwrite('data_set/adapted_meshes.txt',adapted_meshes,'delimiter','\t')
% dlmwrite('data_set/initial_deltas.txt',initial_deltas,'delimiter','\t')
% dlmwrite('data_set/final_deltas.txt',final_deltas,'delimiter','\t')
dlmwrite('data_set/solutions.txt',solutions,'delimiter','\t')
% dlmwrite('data_set/gradients.txt',gradients,'delimiter','\t')







