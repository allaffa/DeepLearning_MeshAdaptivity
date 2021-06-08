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

profile = [-2.0971 , -2.0971 , -2.0971 , -2.0971 , -2.0971 , -2.0971 ,...
       -2.0971 , -2.0971 , -2.0971 , -2.0971 , -2.0971 , -2.0971 ,...
       -2.0971 , -2.0971 , -2.0971 , -2.0971 , -2.0971 , -2.0971 ,...
       -2.0971 , -2.0971 , -2.0971 , -2.0971 , -2.0971 , -2.0971 ,...
       -2.0971 , -2.0971 , -2.0971 , -2.0971 , -2.0971 , -2.0971 ,...
       -2.0971 , -2.0971 , -2.0971 , -2.0971 , -2.0971 , -2.0971 ,...
       -2.0971 , -2.0971 , -2.0971 , -2.0971 , -2.0971 , -2.0971 ,...
       -2.0971 , -2.0971 , -2.0971 , -2.0971 , -2.0971 , -2.0971 ,...
       -2.0971 , -2.0971 , -2.0971 , -2.0971 , -2.0971 , -2.0971 ,...
       -2.0971 , -2.0971 , -2.0971 , -2.0971 , -2.0971 , -2.0971 ,...
       -2.0971 , -2.0971 , -2.0971 , -2.0971 , -2.0971 , -2.0971 ,...
       -2.0971 , -2.0971 , -2.0971 , -2.0971 , -2.0971 , -2.0971 ,...
       -2.0971 , -2.0971 , -2.0971 , -2.0971 , -2.0971 , -2.0971 ,...
       -2.0971 , -2.0971 , -2.0971 , -2.0971 , -2.0971 , -2.0971 ,...
       -2.0971 , -2.0971 , -2.0971 , -2.0971 , -2.0971 , -2.0971 ,...
       -2.0971 , -2.0971 , -2.0971 , -2.0971 , -2.0971 , -2.0971 ,...
       -2.0971 , -2.0971 , -2.0971 , -2.0971 , -2.0971 , -2.0971 ,...
       -2.0971 , -2.0971 , -2.0971 , -2.0971 , -2.0971 , -2.0971 ,...
       -2.0971 , -2.0971 , -2.0971 , -2.0971 , -2.0971 , -2.0971 ,...
       -2.0971 , -2.0955 , -2.0735 , -1.9391 , -1.4462 , -0.25117,...
        1.7548 ,  4.1345 ,  6.1405 ,  7.3355 ,  7.8284 ,  7.9628 ,...
        7.9848 ,  7.9864 ,  7.9864 ,  7.9864 ,  7.9864 ,  7.9864 ,...
        7.9864 ,  7.9864 ,  7.9864 ,  7.9864 ,  7.9864 ,  7.9864 ,...
        7.9864 ,  7.9864 ,  7.9864 ,  7.9864 ,  7.9864 ,  7.9864 ,...
        7.9864 ,  7.9864 ,  7.9864 ,  7.9864 ,  7.9864 ,  7.9864 ,...
        7.9864 ,  7.9864 ,  7.9864 ,  7.9864 ,  7.9864 ,  7.9864 ,...
        7.9864 ,  7.9864 ,  7.9864 ,  7.9864 ,  7.9864 ,  7.9864 ,...
        7.9864 ,  7.9864 ,  7.9864 ,  7.9864 ,  7.9864 ,  7.9864 ,...
        7.9864 ,  7.9864 ,  7.9864 ,  7.9864 ,  7.9864 ,  7.9864 ,...
        7.9864 ,  7.9864 ,  7.9864 ,  7.9864 ,  7.9864 ,  7.9864 ,...
        7.9864 ,  7.9864 ,  7.9864 ,  7.9864 ,  7.9864 ,  7.9864 ,...
        7.9864 ,  7.9864 ,  7.9864 ,  7.9864 ,  7.9864 ,  7.9864 ,...
        7.9864 ,  7.9864 ,  7.9864 ,  7.9864 ,  7.9864 ,  7.9864 ,...
        7.9864 ,  7.9864 ,  7.9864 ]';
  
uniform_grid = linspace(0,1,201)';
[old_physical_mesh] = adaptive_mesh_moving1D_optimization(201, uniform_grid, profile, 0);  
old_deep_mesh = double(python_module.evaluate_model_load(profile)); 

figure()
plot(max(0,uniform_grid), max(0,old_physical_mesh), '-*', 'linewidth', 10);
set(gca, 'fontsize', 45);
xlabel('Uniform mesh coordinates')
ylabel('Adapted mesh coordinates')
axis equal
xlim([0 1])
ylim([0 1])
hold on
plot(max(0,uniform_grid), max(0,old_deep_mesh), '-*', 'linewidth', 10);
legend('Standard adaptivity', 'Deep learning adaptivity');
    
new_profile = profile;
new_profile(30:85) = 7.5;
new_profile(1:90) = smoothdata(new_profile(1:90), 'gaussian', 15);

figure()
plot(max(0,uniform_grid), profile, '-*', 'linewidth', 10);
set(gca, 'fontsize', 45);
ylim([min(profile),max(profile)])
xlabel('Uniform mesh coordinates')
ylabel('Values of the profile')
hold on
plot(max(0,uniform_grid), new_profile, '-*', 'linewidth', 10);
legend('Old shock profile', 'New shock profile');

[new_physical_mesh] = adaptive_mesh_moving1D_optimization(201, uniform_grid, new_profile, 0);
new_deep_mesh = double(python_module.evaluate_model_load(new_profile)); 
new_deep_mesh = sort(new_deep_mesh); 

figure()
plot(max(0,uniform_grid), new_physical_mesh, '-*', 'linewidth', 10);
set(gca, 'fontsize', 45);
xlabel('Uniform mesh coordinates')
ylabel('Adapted mesh coordinates')
axis equal
xlim([0 1])
ylim([0 1])
hold on
plot(max(0,uniform_grid), old_deep_mesh, '-*', 'linewidth', 10);
plot(max(0,uniform_grid), new_deep_mesh, '-*', 'linewidth', 10);
legend('Standard adaptivity', 'OLD - Deep learning adaptivity', 'NEW - Deep learning adaptivity');

new_profile2 = profile;
new_profile2(30:85) = 7.5;
new_profile2(170:185) = 3.2;
new_profile(1:90) = smoothdata(new_profile(1:90), 'gaussian', 15);
new_profile2(170:end) = smoothdata(new_profile2(170:end), 'gaussian', 15);

[new_physical_mesh2] = adaptive_mesh_moving1D_optimization(201, uniform_grid, new_profile2, 0);
new_deep_mesh2 = double(python_module.evaluate_model_load(new_profile2)); 
new_deep_mesh2 = sort(new_deep_mesh2); 

figure()
plot(max(0,uniform_grid), new_physical_mesh2, '-*', 'linewidth', 10);
set(gca, 'fontsize', 45);
xlabel('Uniform mesh coordinates')
ylabel('Adapted mesh coordinates')
axis equal
xlim([0 1])
ylim([0 1])
hold on
plot(max(0,uniform_grid), old_deep_mesh, '-*', 'linewidth', 10);
plot(max(0,uniform_grid), new_deep_mesh2, '-*', 'linewidth', 10);
legend('Standard adaptivity', 'OLD - Deep learning adaptivity', 'NEW - Deep learning adaptivity');

