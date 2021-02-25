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
niter_mmpde  = 4;
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

%Initialize physical mesh as uniform
physical_mesh = linspace( xa, xb, nnodes )';

num_problems = 500000;
range_value = [-10; 10];
max_num_shocks = 4;

index_counter = 1;
shock_coordinate = xa;

shocks = zeros(num_problems, nnodes); 
final_meshes = zeros(num_problems, nnodes); 

value = range_value(1) + ( range_value(2)-range_value(1) )*rand(1,1);

for problem_index = 1:num_problems

    num_shocks = randi(max_num_shocks,1,1); 
    profile = zeros(size(comp_mesh));
    
    shocks_position = xa+(xb-xa)*rand(1,num_shocks);
    shocks_position = sort(shocks_position);
    
    shock_counter = 1;
    
    for index_coordinate =1:nnodes
        if( (shock_counter > num_shocks) || (comp_mesh(index_coordinate) <= shocks_position(shock_counter)) )
            profile(index_coordinate) = value;
        else
            shock_counter = shock_counter + 1;
            value = range_value(1) + ( range_value(2)-range_value(1) )*rand(1,1);
            profile(index_coordinate) = value;
        end
    end
    
    num_profile_smoothing = rand(1,3);
    
    for iter = 0:num_profile_smoothing
        profile = smoothdata(profile, 'gaussian', ceil(ksmooth*length(comp_mesh)));
    end
    
    [physical_mesh,u0,u1,omega] = adapt_mesh(nnodes,xa,xb,comp_mesh,profile,profile,physical_mesh,niter_mmpde,eps_omega,avrg,smoothing,nsmooth,ksmooth,interp_method);
    grid = physical_mesh';
    
    index_counter = index_counter + 1;
    
    shocks(problem_index,:) = profile';
    final_meshes(problem_index,:) = grid;    

end


dlmwrite('shocks.txt', shocks,'delimiter','\t')
dlmwrite('final_meshes.txt',final_meshes,'delimiter','\t')




        
