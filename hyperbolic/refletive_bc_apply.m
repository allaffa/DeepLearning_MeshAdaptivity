function q = reflet_bc_apply(q,nghost)

nvar = size(q,1);
ncell = size(q,2);

q(1,1:1:nghost) =  q(1,nghost+nghost:-1:1+nghost);
q(2,1:1:nghost) = -q(2,nghost+nghost:-1:1+nghost);
q(3,1:1:nghost) =  q(3,nghost+nghost:-1:1+nghost);

q(1,ncell-nghost+1:ncell) =   q(1,ncell-nghost:-1:ncell-2*nghost+1);
q(2,ncell-nghost+1:ncell) =  -q(2,ncell-nghost:-1:ncell-2*nghost+1);
q(3,ncell-nghost+1:ncell) =   q(3,ncell-nghost:-1:ncell-2*nghost+1);
end



