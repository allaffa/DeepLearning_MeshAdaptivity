function [data] = exact_sedov_1d(gam, g_x)
%here r is the designated point, i.e., cell center c_x

c_x = ( g_x(1:end-1) + g_x(2:end) )/2;

% Sedov uniform distributed case
w=0;

%initial values
j=1;

a=(j+2-w)*(gam+1)/4.0;
b=(gam+1)/(gam-1);
c=(j+2-w)*gam/2.0;
d=(j+2-w)*(gam+1)/( (j+2-w)*(gam+1) - 2.0*(2.0+j*(gam-1.0)) );
e=(2+j*(gam-1.0))/2.0;

%%% checkthis alp0
alp0 = 2.0/(j+2-w);
% alp0 = 2.0/(j-2)
alp2 = -(gam-1.0)/( 2.0*(gam-1.0) + j - gam*w  );
alp1 = (j+2.0-w)*gam/(2.0+j*(gam-1.0)) * ( 2.0*(j*(2.0-gam)-w)/gam/(j+2.0-w)^2 - alp2 );
%%% checkthis alp3
alp3 = (j-w) / ( 2.0*(gam-1.0) + j -gam*w);
alp4 = (j+2-w)*(j-w)/(j*(2.0-gam)-w)*alp1;
alp5 = ( w*(1.0+gam)-2.0*j )/ ( j*(2.0-gam) -w );

% desired values
r2=0.5;
rho2=6.0;
v2=0.277778;
E2=0.0385802;
p2=0.0925926;

n=length(c_x);
v=zeros(n,1);
p=zeros(n,1);
rho=zeros(n,1);
sie=zeros(n,1);

% compute v,rho,p,e after compute V for given r(i) location
for i=1:200
  func = @(x)(r2*(a*x).^(-alp0).*(b*(c*x-1.0)).^(-alp2).*(d*(1-e*x)).^(-alp1) -c_x(i));
  V = fzero(  func, 0.517599999999 );
%   V = fzero(  func, 0.5175999999999999 );
% V = fzero( @(x)( r2*(a*x)^(-alp0)*(b*(c*x-1.0))^(-alp2) -r(i) ), 2.0)
 lambda = c_x(i)/r2;
 
 x1 = a*V;
 x2 = b*(c*V-1.0);
 x3 = d*(1-e*V);
 x4 = b*(1-c/gam*V);
 
 v(i) = v2*x1*lambda;
 rho(i) = rho2*x1^(alp0*w)*x2^(alp3+alp2*w)*x3^(alp4+alp1*w)*x4^alp5;
 p(i) = p2*x1^(alp0*j)*x3^(alp4+alp1*(w-2))*x4^(1+alp5);
 sie(i) = p(i)/(gam-1.0)/rho(i);
end

% % interpolate with pchip
% v(1:11,1) = pchip(r(12:200,1), v(12:200,1), r(1:11,1) );
rho(1:17,1) = pchip(c_x(18:200,1), rho(18:200,1), c_x(1:17,1) );
% p(1:11,1) = pchip(r(18:200,1), p(18:200,1), r(1:11,1) );
% sie(1:11,1) = pchip(r(18:200,1), sie(18:200,1), r(1:11,1) );

% extrapolate
v(1:17,1) =interp1(c_x(18:200,1), v(18:200,1), c_x(1:17,1),'pchip','extrap' );
% rho(1:17,1) = interp1(r(18:200,1), rho(18:200,1), r(1:17,1),'linear','extrap' );
p(1:17,1) = interp1(c_x(18:200,1), p(18:200,1), c_x(1:17,1),'pchip','extrap' );
sie(1:17,1) = interp1(c_x(18:200,1), sie(18:200,1), c_x(1:17,1),'pchip','extrap' );

% interpolate grid vel from cell to grid
v = pchip(c_x(1:200,1), v(1:200,1), g_x(1:201));
v(1) = 0.0;


% initial values outside of shock
rho(201:240) = 1.0;
p(201:240) = 0.0;
sie(201:240) = 0.0;
v(201:241) = 0.0;

data.x = g_x;
data.c_x = c_x;
data.rho = rho;
data.u = v;
data.e = sie;
data.P = p;
