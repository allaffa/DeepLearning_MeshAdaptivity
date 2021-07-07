%Physical domain boundaries
xa = 0.0;
xb = 1.0;

%Parameters
smooth_opt   = ["none","weightavrg","gaussian"];
avrg_opt     = ["arithmetic","harmonic","geometric"];
interp_opt   = ["linear","pchip"];
nnodes       = 201;
time_start   = 0.0;
velocity     = 1.e-2;
diffusion    = 0.e-5;
niter_mmpde  = 1;
eps_nonlin   = 1e-5;
niter_coupl  = 2;
smoothing    = smooth_opt(3);
nsmooth      = 2;
ksmooth      = 1.e-1;
avrg         = avrg_opt(1); 
numtimesteps = 16;
timestep     = 0.5;
time 	       = time_start;
eps_omega    = 1.e-2;
interp_method= interp_opt(1);

%Computational mesh is always uniform; ksi e [0,1]
comp_mesh = linspace( 0, 1, nnodes)';

num_problems = 1000;
range_value = [-10; 10];
max_num_shocks = 4;

index_counter = 1;
shock_coordinate = xa;

shocks = zeros(num_problems, nnodes); 
final_meshes = zeros(num_problems, nnodes); 

value = range_value(1) + ( range_value(2)-range_value(1) )*rand(1,1);

num_shocks = randi(max_num_shocks,1,1); 

shocks_position = xa+(xb-xa)*rand(1,num_shocks);
shocks_position = sort(shocks_position);
    

for problem_index = 1:num_problems
    
    %Initialize physical mesh as uniform
    physical_mesh = linspace( xa, xb, nnodes )';

    profile_function = @(x, t, v) 0.0 + 0.6 * ( ( (x >= 0.3 + v * t) & (x <= 0.45 + v * t ) ) ) + 0.4 * ( ( (x >= 0.35 + v * t) & (x <= 0.65 + v * t ) ) );
    
    profile = profile_function(physical_mesh, 0, 0);
    
    num_profile_smoothing = 1;
    
    for iter = 0:num_profile_smoothing
        profile = smoothdata(profile, 'gaussian', ceil(ksmooth*length(physical_mesh)));
    end
    
    [physical_mesh,u0,u1,omega] = adapt_mesh(nnodes,xa,xb,comp_mesh,profile,profile,physical_mesh,physical_mesh,niter_mmpde,eps_nonlin,eps_omega,avrg,smoothing,nsmooth,ksmooth,interp_method);
                                 
    grid = physical_mesh';
    
    index_counter = index_counter + 1;
    
    shocks(problem_index,:) = profile';
    final_meshes(problem_index,:) = grid;    

end


dlmwrite('shocks.txt', shocks,'delimiter','\t')
dlmwrite('final_meshes.txt',final_meshes,'delimiter','\t')




        
