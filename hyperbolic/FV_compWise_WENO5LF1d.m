function res = FV_compWise_WENO5LF1d(gamma,lambda,dt,q,dx)

nx=size(q,2); R=3; I=R:(nx-R); % R: stencil size 
LF=zeros(size(q)); res=zeros(size(q));

%% Right Flux
% Choose the positive fluxes, 'v', to compute the left cell boundary flux:
% $u_{i+1/2}^{-}$
vmm = q(:,I-2);
vm  = q(:,I-1);
v   = q(:, I );
vp  = q(:,I+1);
vpp = q(:,I+2);

% Polynomials
p0n = (2*vmm - 7*vm + 11*v)/6;
p1n = ( -vm  + 5*v  + 2*vp)/6;
p2n = (2*v   + 5*vp - vpp )/6;

% Smooth Indicators (Beta factors)
B0n = 13/12*(vmm-2*vm+v  ).^2 + 1/4*(vmm-4*vm+3*v).^2; 
B1n = 13/12*(vm -2*v +vp ).^2 + 1/4*(vm-vp).^2;
B2n = 13/12*(v  -2*vp+vpp).^2 + 1/4*(3*v-4*vp+vpp).^2;

% Constants
d0n = 1/10; d1n = 6/10; d2n = 3/10; epsilon = 1e-6;

% Alpha weights 
alpha0n = d0n./(epsilon + B0n).^2;
alpha1n = d1n./(epsilon + B1n).^2;
alpha2n = d2n./(epsilon + B2n).^2;
alphasumn = alpha0n + alpha1n + alpha2n;

% ENO stencils weigths
w0n = alpha0n./alphasumn;
w1n = alpha1n./alphasumn;
w2n = alpha2n./alphasumn;

% Numerical Flux at cell boundary, $u_{i+1/2}^{-}$;
hn = w0n.*p0n + w1n.*p1n + w2n.*p2n;

%% Left Flux 
% Choose the negative fluxes, 'u', to compute the left cell boundary flux:
% $u_{i-1/2}^{+}$ 
umm = q(:,I-1);
um  = q(:, I );
u   = q(:,I+1);
up  = q(:,I+2);
upp = q(:,I+3);

% Polynomials
p0p = ( -umm + 5*um + 2*u  )/6;
p1p = ( 2*um + 5*u  - up   )/6;
p2p = (11*u  - 7*up + 2*upp)/6;

% Smooth Indicators (Beta factors)
B0p = 13/12*(umm-2*um+u  ).^2 + 1/4*(umm-4*um+3*u).^2; 
B1p = 13/12*(um -2*u +up ).^2 + 1/4*(um-up).^2;
B2p = 13/12*(u  -2*up+upp).^2 + 1/4*(3*u -4*up+upp).^2;

% Constants
d0p = 3/10; d1p = 6/10; d2p = 1/10; epsilon = 1e-6;

% Alpha weights 
alpha0p = d0p./(epsilon + B0p).^2;
alpha1p = d1p./(epsilon + B1p).^2;
alpha2p = d2p./(epsilon + B2p).^2;
alphasump = alpha0p + alpha1p + alpha2p;

% ENO stencils weigths
w0p = alpha0p./alphasump;
w1p = alpha1p./alphasump;
w2p = alpha2p./alphasump;

% Numerical Flux at cell boundary, $u_{i-1/2}^{+}$;
hp = w0p.*p0p + w1p.*p1p + w2p.*p2p;

%% Compute finite volume residual term, df/dx.
LF(:,I) = 0.5*(F(hn, gamma)+F(hp, gamma)-abs(lambda).*(hp-hn)); % Lax friedrichs flux
% for j = I % for all faces of the domain cells
%     res(:, j ) = res(:, j ) + LF(:,j)/dx;
%     res(:,j+1) = res(:,j+1) - LF(:,j)/dx;
% end % or alternatively:
res(:,I) = (LF(:,I)-LF(:,I-1))/dx; % L = -df(q)/dx.

% Flux contribution of the LEFT MOST FACE: left face of cell j=1.
res(:,3) = res(:,3)-LF(:,3)/dx;
 
% Flux contribution of the RIGHT MOST FACE: right face of cell j=nx-1.
res(:,nx-2)=res(:,nx-2)+LF(:,nx-2)/dx;
end

% Compute flux vector
function flux = F(q, gamma)
    
    % primary properties
    rho=q(1,:); u=q(2,:)./rho; E=q(3,:); p=(gamma-1)*(E-0.5*rho.*u.^2);
    
    % flux vector of conserved properties
    flux=[rho.*u; rho.*u.^2+p; u.*(E+p)];
end