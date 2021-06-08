function [res] = hyperbolic_PDE_solver_WENO5_nonuniform(x, sol, c, dt)

    nx=size(sol,1); R=3; I=R:(nx-R); % R: stencil size 
    LF=zeros(size(sol)); res=zeros(size(sol));
    DX = zeros(size(x));
    DX(1) = x(2)-x(1);
    DX(end) = x(end) - x(end-1);
    
    for i =2:length(x)-1
        DX(i) = 0.5 * ((x(i)-x(i-1))+(x(i+1)-x(i)));
    end

    %% Right Flux
    % Choose the positive fluxes, 'v', to compute the left cell boundary flux:
    % $u_{i+1/2}^{-}$
    vmm = sol(I-2);
    vm  = sol(I-1);
    v   = sol(I );
    vp  = sol(I+1);
    vpp = sol(I+2);
    
    dxmm = DX(I-2);
    dxm = DX(I-1);
    dx = DX(I);
    dxp = DX(I+1);
    dxpp = DX(I+2);

    % Polynomials
    
    A = 2/3 * dxm./(dxmm+dxm);
    B = -1/3 * (1 - dxm./(dxm+dx)) - 4/3 * dxm./(dxmm+dxm) - 2/3 * dxmm./(dxmm+dxm);
    C = 2 - 1/3* dxm./(dxm+dx);
    p0n = A .* vmm + B .* vm + C .* v;
    
    D = -1/3 * (1-dxm./(dxm+dx));
    E = -1/3 * dxm./(dxm+dx) + 2/3 + 2/3 * (1-dx./(dx+dxp));
    F = 2/3 * dx./(dx+dxp);
    p1n = D .* vm  + E .* v  + F .* vp; 
    
    G = 2/3 * (1 - dx./(dx+dxp));
    H = 2/3 * dx./(dx+dxp) + 2/3 - 1/3 * (1 - dxp./(dxp+dxpp));
    J = -1/3 * dxp./(dxp + dxpp);    
    p2n =  G .* v + H .* vp + J .* vpp;
    
    % Smooth Indicators (Beta factors)
    
    A1 = 2 * dxm./(dxmm + dxm);
    B1 = 2 * (1-dxm./(dxm+dx)) - 4 * dxm./(dxmm+dxm) - 2 * dxmm./(dxmm+dxm);
    C1 = 2 * dxm./(dxm+dx);
    
    A2 = 2 * dxm ./(dxmm+dxm);
    B2 = -2 * (1-dxm./(dxm+dx)) - 4 * dxm./(dxmm+dxm) -2 * dxmm./(dxmm+dxm);
    C2 = 4 - 2 * dxm./(dxm+dx);
    
    B0n = 13/12*(A1.*vmm + B1.*vm+C1.*v  ).^2 + 1/4*( A2.*vmm + B2.*vm + C2.*v).^2; 
      
    D1 = 2*(1-dxm./(dxm+dx));
    E1 = 2*dxm./(dxm+dx) -4 + 2*(1-dx./(dx+dxp));
    F1 = 2*dx./(dx+dxp);
    
    D2 = -2 * (1-dxm./(dxm+dx)) + 2 * (1-dx./(dx+dxp));
    E2 = 0.0;
    F2 = -2 * dxm./(dxm+dx) + 2 * dx./(dx+dxp);
    
    B1n = 13/12*(D1.*vm + E1.*v + F1.*vp ).^2 + 1/4*(D2.*vm + E2.* v + F2.*vp).^2;
    
    G1 = 2*(1 - dx./(dx+dxp));
    H1 = 2*(dx./(dx+dxp)) - 4 + 2*(1 - dxp./(dxp+dxpp));
    J1 = 2 * (dxp./(dxp+dxpp));
    
    G2 = -6 * (1 - dx./(dx+dxp));
    H2 = -6 * dx./(dx+dxp) + 8 -2 * (1 - dxp./(dxp+dxpp));
    J2 = -2 * dxp./(dxp+dxpp);
    
    B2n = 13/12*(G1.*v + H1.*vp + J1.*+vpp).^2 + 1/4*(G2.*v +H2.*vp + J2.*vpp).^2;    

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
    umm = sol(I-1);
    um  = sol(I);
    u   = sol(I+1);
    up  = sol(I+2);
    upp = sol(I+3);
    
    dxmm = DX(I-1);
    dxm = DX(I);
    dx = DX(I+1);
    dxp = DX(I+2);
    dxpp = DX(I+3); 
    
    % Polynomials
    
    A = -1/3 * (1-dxmm./(dxmm+dxm));
    B = -1/3 * dxm./(dxm+dx) + 2/3 + 2/3 * (1-dx./(dx+dxp));
    C = 2/3 * dx./(dx+dxp);
    p0p = A .* umm + B .* um + C .* u ;    
    
    D = 2/3 * (1 - dxm./(dxm+dx));
    E = ( 2/3 * dxm./(dxm+dx) + 2/3 - 1/3 * (1 - dx./(dx+dxp)) );
    F = -1/3 * dx./(dx + dxp);    
    p1p = D .* um + E .* u  + F .* up;
    
    G = 2 - 1/3 * dx./(dx+dxp);
    H = -1/3 * (1 - dxp./(dxp+dxpp)) - 4/3 * dxp./(dx+dxp) - 2/3 * dx./(dx+dxp);
    J = 2/3* dxp./(dxp+dxpp);
    p2p = G .* u + H .* up + J .* upp;
    
    % Smooth Indicators (Beta factors)
    
    A1 = 2 * dxm./(dxmm + dxm);
    B1 = 2 * (1-dxm./(dxm+dx)) - 4 * dxm./(dxmm+dxm) - 2 * dxmm./(dxmm+dxm);
    C1 = 2 * dxm./(dxm+dx);
    
    A2 = 2 * dxm ./(dxmm+dxm);
    B2 = -2 * (1-dxm./(dxm+dx)) - 4 * dxm./(dxmm+dxm) -2 * dxmm./(dxmm+dxm);
    C2 = 4 - 2 * dxm./(dxm+dx);
    
    B0p = 13/12*(A1.*umm + B1.*um+C1.*u  ).^2 + 1/4*( A2.*umm + B2.*um + C2.*u).^2; 
      
    D1 = 2*(1-dxm./(dxm+dx));
    E1 = 2*dxm./(dxm+dx) -4 + 2*(1-dx./(dx+dxp));
    F1 = 2*dx./(dx+dxp);
    
    D2 = -2 * (1-dxm./(dxm+dx)) + 2 * (1-dx./(dx+dxp));
    E2 = 0.0;
    F2 = -2 * dxm./(dxm+dx) + 2 * dx./(dx+dxp);
    
    B1p = 13/12*(D1.*um + E1.*u + F1.*up ).^2 + 1/4*(D2.*um + E2.* u + F2.*up).^2;
    
    G1 = 2*(1 - dx./(dx+dxp));
    H1 = 2*(dx./(dx+dxp)) - 4 + 2*(1 - dxp./(dxp+dxpp));
    J1 = 2 * (dxp./(dxp+dxpp));
    
    G2 = -6 * (1 - dx./(dx+dxp));
    H2 = -6 * dx./(dx+dxp) + 8 -2 * (1 - dxp./(dxp+dxpp));
    J2 = -2 * dxp./(dxp+dxpp);
    
    B2p = 13/12*(G1.*u + H1.*up + J1.*+upp).^2 + 1/4*(G2.*u +H2.*up + J2.*upp).^2;    
    
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
    LF(I) = 0.5*(c*hn+c*hp-abs(dt).*(hp-hn)); % Lax friedrichs flux
    % for j = I % for all faces of the domain cells
    %     res(:, j ) = res(:, j ) + LF(:,j)/dx;
    %     res(:,j+1) = res(:,j+1) - LF(:,j)/dx;
    % end % or alternatively:
    res(I) = (LF(I)-LF(I-1))./DX(I); % L = -df(sol)/dx.

    % Flux contribution of the LEFT MOST FACE: left face of cell j=1.
    res(3) = res(3)-LF(3)./DX(3);

    % Flux contribution of the RIGHT MOST FACE: right face of cell j=nx-1.
    res(nx-2)=res(nx-2)+LF(nx-2)./(dx(end-1));

end