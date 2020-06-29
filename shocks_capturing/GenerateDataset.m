n0 = 41;
x0 = 0; 
x1 = 1;
reference_grid = linspace(x0, x1, n0)'; 
num_problems = 50000;
range_value = [-10; 10];
max_num_shocks = 3;

index_counter = 1;
shock_coordinate = x0;
kernel_width = 5;

shocks = zeros(num_problems, n0); 
final_meshes = zeros(num_problems, n0); 

value = range_value(1) + ( range_value(2)-range_value(1) )*rand(1,1);

for problem_index = 1:num_problems

    num_shocks = randi(max_num_shocks,1,1); 
    profile = zeros(size(reference_grid));
    
    shocks_position = x0+(x1-x0)*rand(1,num_shocks);
    shocks_position = sort(shocks_position);
    
    shock_counter = 1;
    
    for index_coordinate =1:n0
        if( (shock_counter > num_shocks) || (reference_grid(index_coordinate) <= shocks_position(shock_counter)) )
            profile(index_coordinate) = value;
        else
            shock_counter = shock_counter + 1;
            value = range_value(1) + ( range_value(2)-range_value(1) )*rand(1,1);
            profile(index_coordinate) = value;
        end
    end
    
    [adapted_mesh] = adaptive_mesh_moving1D_optimization(reference_grid, profile, kernel_width);
    grid = adapted_mesh';
    
    index_counter = index_counter + 1;
    
    shocks(problem_index,:) = profile';
    final_meshes(problem_index,:) = grid;    

end


dlmwrite('shocks.txt', shocks,'delimiter','\t')
dlmwrite('final_meshes.txt',final_meshes,'delimiter','\t')




        
