function [ u, x, y, t, P ] ...
    = solve_st_neumann( u0, T, area, h, k, n )
%SOLVE_ST_NEUMANN Solves the diffusion equation inside a rectangle in R^2
% with zero Neumann boundary conditions using Brownian motion.
%   Args:
%     u0        struct
%                   .fun    Function returning initial values.
%                   .supp   Support of the function.
%                   .m      (opional) Maximal value of the function.
%               or a function handle corresponding to u0.fun above.
%     T         Length of time interval.
%     area      2x2 matrix representing rectangle Q as a cartesian
%                 product of its rows treated as closed intervals.
%     h         Spatial step size, represented as a 2-element vector or
%                 a scalar.
%     k         Temporal step size.
%     n         Number of samples from initial distribution.
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

itg = integral2(u0, supp(1,1), supp(1,2), supp(2,1), supp(2,2));

if hasm
    P0 = sample_rejection(u0,n,supp,m);
else
    P0 = sample_rejection(u0,n,supp);
end
[P, t] = diffuse_neumann(P0,area,T,k);

[u_tmp, x, y] = normalize(P(:,:,1), itg, n, area, h);
lt = length(t); lx = length(x); ly = length(y);
u = zeros(lx, ly, lt);
u(:,:,1) = u_tmp;
for i=2:lt
    u(:,:,i) = normalize(P(:,:,i), itg, n, area, h);
end

end
