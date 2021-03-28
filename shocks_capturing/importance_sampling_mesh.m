function mesh = importance_sampling_mesh(initial_mesh, solution, gamma_val)

    solution = smoothdata(solution, 'gaussian', ceil(0.1*length(solution)));

    final_mesh = zeros(size(initial_mesh));
    n_iters = 1;
    
    mesh_size = size(initial_mesh);
    standard_deviation =  zeros(mesh_size);
    
    spacing = diff(initial_mesh);

    for entry=2:mesh_size-1
        standard_deviation(entry) = min([spacing(entry-1) spacing(entry)]);
    end
    
    for iter = 1:n_iters
        
        state = initial_mesh;

        no_action = zeros(mesh_size);
        [~, prev_reward]= compute_mesh_uniform_gradient_norm(state, solution, no_action);
        
        for step_count = 1:5
            
            spacing = diff(initial_mesh);
            for entry=2:mesh_size-1
                standard_deviation(entry) = min([spacing(entry-1) spacing(entry)])/2;
            end            
            
            reward = prev_reward-1;
            count = 0;
            
            while reward < prev_reward && count < 100
                action = zeros(mesh_size);
                action(1) = 0.0;
                for entry=2:mesh_size-1
                    action(entry) = standard_deviation(entry) * ( rand - 0.5 );
                end
                action(end) = 0.0;
                [next_state, reward]= compute_mesh_uniform_gradient_norm(state, solution, action);
%                 display(strcat("prev_reward: ", num2str(prev_reward), " reward: ", num2str(reward)))
                count = count + 1;
            end
            
            if reward > prev_reward
                prev_reward = reward;
                state = next_state;
            end
            
        end
       
        final_mesh = final_mesh + state;
        
    end
        
    mesh = final_mesh/n_iters;
        
end



function [current_state, reward] = compute_mesh_uniform_gradient_norm(initial_state, solution, action)

    current_state = initial_state + action;
    current_state = sort(current_state);
    current_solution = interp1(initial_state,solution,current_state,"linear",'extrap');
   
    grad_norm = abs(diff(current_solution));
    mesh_degeneracy = diff(initial_state);
   
    reward = std((grad_norm.*mesh_degeneracy));
    reward = -log(reward);
    
end



