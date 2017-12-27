function [ u, x, y, t, P ] = solve_st( u0, bdr, T, area, h, k, n, boundary)
%SOLVE_ST Solves the diffusion equation inside a region in R^2
% with zero Dirichlet/Neumann boundary conditions using Brownian motion.
%   Args:
%     u0        struct
%                   .fun    Function returning initial values or []
%                           for delta_0 distribution.
%                   .supp   Support of the function.
%                   .m      (opional) Maximal value of the function.
%               or a function handle corresponding to u0.fun above.
%     bdr       Boolean function returning 1 if a point lies inside
%                 the region.
%     T         Length of time interval.
%     area      2x2 matrix representing solution bounding box
%                 as a cartesian product of its rows treated as closed
%                 intervals.
%     h         Spatial step size, represented as a 2-element vector or
%                 a scalar.
%     k         Temporal step size.
%     n         Number of samples from initial distribution.
%     boundary  Name of the boundary condition. Available names:
%                 'None'        No boundary. Let the region be whole R^2.
%                 'Dirichlet'   Absorbing boundary. The region is described
%                                 by the function `bdr`.
%                 'Neumann'     Reflecting boundary. The region is a
%                                 rectangle bounded by `area`.
%
%   Returns:
%     u         Tensor of approximate solutions. First two coordinates are
%                 spatial, while the last is temporal.
%     x,y       Coordinates of spatial mesh nodes.
%     t         Mesh of time discretization.
%     P         Trajectories of each particle as 2xnxm tensor, where
%                 m is the length of t.

hasm = 0;
if isstruct(u0)
    if isfield(u0,'m')
        m = u0.m;
        hasm = 1;
    end
    supp = u0.supp;
    u0 = u0.fun;
else
    supp = area;
end

d = size(area,1);

if isempty(u0)
    P0 = zeros(d,n);
    itg = 1;
else
    if hasm
        P0 = sample_rejection(u0,n,supp,m);
    else
        P0 = sample_rejection(u0,n,supp);
    end
    itg = integral2(u0, supp(1,1), supp(1,2), supp(2,1), supp(2,2));
end

if strcmp(boundary, 'None')
    [P, t] = diffuse_nobd(P0,T,k);
elseif strcmp(boundary, 'Dirichlet')
    [P, t] = diffuse_dirichlet(P0,bdr,T,k);
    P2 = P;
    b = squeeze(bdr(P));
    for j=1:length(t)
        P2(:,~b(:,j),j) = NaN;
    end
elseif strcmp(boundary, 'Neumann')
    [P, t] = diffuse_neumann(P0,area,T,k);
else
    error('Error: boundary condition %s not implemented', boundary);
end

[~, x, y] = normalize(P(:,:,1), itg, n, area, h);
lt = length(t); lx = length(x); ly = length(y);
u = zeros(lx, ly, lt);
if strcmp(boundary, 'Dirichlet')
    for i=1:lt
        u(:,:,i) = normalize(P2(:,:,i), itg, n, area, h);
    end
else
    for i=1:lt
        u(:,:,i) = normalize(P(:,:,i), itg, n, area, h);
    end
end

end

