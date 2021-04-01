%Physical domain boundaries
xa = 0.0;
xb = 1.0;

%Parameters
smooth_opt      = ["none","weightavrg","gaussian"];
avrg_opt        = ["arithmetic","harmonic","geometric"];
interp_opt      = ["linear","pchip"];
nnodes          = 201;
time_start      = 0.0;
velocity        = 1.e-2;
diffusion       = 0.e-5;
niter_mmpde     = 128;
eps_nonlin      = 1e-5;
niter_coupl     = 1;
smoothing       = smooth_opt(2);
nsmooth_omega   = 4;
nsmooth_profile = 0;
ksmooth_omega   = 1.e-1;
ksmooth_profile = 1.e-1;
avrg            = avrg_opt(1); 
numtimesteps    = 16;
timestep        = 0.5;
time 	        = time_start;
eps_omega       = 1.e-2;
interp_method   = interp_opt(1);

%Computational mesh is always uniform; ksi e [0,1]
comp_mesh = linspace( 0, 1, nnodes)';

%Initialize physical mesh as uniform
physical_mesh = linspace( xa, xb, nnodes )';

num_problems = 1;
range_value = [0; 1];
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

    %%Smooth profile
    if smoothing=="weightavrg" profile = smooth_fun(profile,nsmooth_profile,1);, end
    if smoothing=="gaussian"   profile = smoothdata(profile, 'gaussian', ceil(ksmooth_profile*length(physical_mesh)));,end    

    %MMPDE solver
    [physical_mesh,u0,u1,omega] = adapt_mesh(nnodes,xa,xb,comp_mesh,profile,profile,physical_mesh,physical_mesh,niter_mmpde,eps_nonlin,eps_omega,avrg,smoothing,nsmooth_omega,ksmooth_omega,interp_method);
                                 
    grid = physical_mesh';
    
    index_counter = index_counter + 1;
    
    shocks(problem_index,:) = profile';
    final_meshes(problem_index,:) = grid;    

    %CV dx
    x = physical_mesh;
    delx = zeros(nnodes,1);
    delx(1,1) = x(2)-x(1);
    for i = 2:nnodes-1
      delx(i,1) = 0.5*(x(i+1)-x(i-1));
    end
    delx(nnodes,1) = x(nnodes)-x(nnodes-1);
    delx = delx/max(delx);

    %figure()
    %plot(physical_mesh,profile,'-','LineWidth',3)
    %hold on
    %plot(physical_mesh,delx,'-')
    %hold on
    %plot(physical_mesh,omega,'--')
    %hold on
    %plot(comp_mesh,physical_mesh,'.')
    %%plot(physical_mesh,0.5,'s')
    %%xlim([xa,xb]);
    %%ylim([-0.1,1.1]);
    %legend('\phi','{\Delta x}','{\omega}','{\xi}');
    

end


dlmwrite('shocks.txt', shocks,'delimiter','\t')
dlmwrite('final_meshes.txt',final_meshes,'delimiter','\t')




        
