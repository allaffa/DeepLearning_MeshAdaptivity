function [flag] = substep_dtnk_new(adapt_mesh, mesh)
n = length( adapt_mesh );

% eq 160 of Adaptive Moving Mesh Method 2011
flag = 1;
for i=2:n-1
 if( 0.5*(mesh(i-1)+mesh(i)) > adapt_mesh(i) || adapt_mesh(i) > 0.5*(mesh(i)+mesh(i+1)) )
   flag = 0;
   break;
 end
end
