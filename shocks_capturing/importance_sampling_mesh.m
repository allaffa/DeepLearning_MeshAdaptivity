function mesh = importance_sampling_mesh(initial_mesh, solution, standard_deviation)

    final_mesh = zeros(size(initial_mesh));
    n_iters = 1000;
    
    mesh_size = size(initial_mesh);
  
    for iter = 1:n_iters
        %print(strcat("Iteration: ",num2str(iter)));
        state = initial_mesh;

        no_action = zeros(mesh_size);
        [~, reference_reward]= compute_mesh_uniform_gradient_norm(state, solution, no_action);
        
        for count_proposal = 1:3
                action = standard_deviation * ( rand(mesh_size) - 0.5 );
                action(1) = 0.0;
                action(end) = 0.0;
                [next_state, reward]= compute_mesh_uniform_gradient_norm(state, solution, action);
                if reward < reference_reward
                    break;
                end
        end

        state = next_state;        
        final_mesh = final_mesh + state;
        
    end
        
    mesh = final_mesh/n_iters;
        
end



function [current_state, reward] = compute_mesh_uniform_gradient_norm(initial_state, solution, action)

    current_state = initial_state + action;
    current_state = sort(current_state);
    current_solution = interp1(initial_state,solution,current_state,"linear",'extrap');
   
    grad_norm = diff(current_solution);
   
    reward = std(grad_norm);
    reward = -log(reward);
    
end



