pyExec	= '/opt/anaconda3/bin/python3';

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

n0 = 201;
x0 = 0; 
x1 = 1;
physical_grid = linspace(x0, x1, n0)'; 
uniform_mesh = physical_grid; 

time_steps = 200;
delta_shift = (x1-x0)/time_steps;

index_counter = 1;
shock_coordinate = x0;
kernel_width = n0;

value = 1.0;

figure
pause

for time_index = 1:time_steps

    solution = zeros(size(physical_grid));
    
    index_coordinate = 1;
    while index_coordinate < n0
        if( physical_grid(index_coordinate) <= shock_coordinate )
            solution(index_coordinate) = value;
        else
            solution(index_coordinate) = 0.0;
        end
        index_coordinate = index_coordinate+1;
    end

    [adapted_mesh] = adaptive_mesh_moving1D_optimization(physical_grid, solution, kernel_width);
    deep_mesh = double(python_module.evaluate_model_load(solution));
    
    subplot(1,2,1)
    plot(physical_grid, solution, '*-', 'LineWidth', 5)
    xlim([x0, x1])
    xlabel('domain')
    ylabel('profile value')
    title('Shock wave')
    set(gca, 'fontsize', 36)
    subplot(1,2,2)
    plot(adapted_mesh, '*-', 'LineWidth', 8)
    hold on
    plot(deep_mesh, 'o-', 'LineWidth', 3)
    plot(uniform_mesh, 'LineWidth', 3)
    xlabel('node index')
    ylabel('node coordinate')
    title('Mesh')
    set(gca, 'fontsize', 36)
    legend('Standard Adaptive', 'Deep Learning', 'Uniform', 'location', 'southeast')
    hold off
    xlim([1 n0])
    ylim([x0 x1])
    pause(0.1)
    
    grid = adapted_mesh';
    
    index_counter = index_counter + 1;
    shock_coordinate = shock_coordinate + delta_shift;

end


