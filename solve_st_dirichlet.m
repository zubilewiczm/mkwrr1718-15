function [ u, x, y, t, P, ep, et ] ...
    = solve_st_dirichlet( u0, bdr, T, area, h, k, n )
%SOLVE_ST_DIRICHLET Solves the diffusion equation inside a region in R^2
% with zero Dirichlet boundary conditions using Brownian motion.
%   Args:
%     u0        struct
%                   .fun    Function returning initial values.
%                   .supp   Support of the function.
%                   .m      (opional) Maximal value of the function.
%               or a function handle corresponding to u0.fun above.
%     bdr       Boolean function returning 1 if a point lies outside
%                 the region.
%     T         Length of time interval.
%     area      2x2 matrix representing solution bounding box
%                 as a cartesian product of its rows treated as closed
%                 intervals.
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
%     ep        Exit points, ie. 2xn matrix, where each column is the
%                 position of the nth particle at the time when it first
%                 leaves the region.
%     et        Exit times, a vector such that t(et) gives the first exit
%                 time for each particle.

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
[P, t] = diffuse_nobd(P0,T,k);

[~, x, y] = normalize(P(:,:,1), itg, n, area, h);
lt = length(t); lx = length(x); ly = length(y);
u = zeros(lx, ly, lt);

where  = squeeze(bdr(P));
oldidx = ones(1,n);
ep     = NaN*zeros(2,n);
et     = NaN*zeros(1,n);
for i=1:lt
    idx = all(where(:,1:i), 2).';
    u(:,:,i) = normalize(P(:,idx,i), itg, n, area, h);
    ep(:,oldidx & ~idx) = P(:,oldidx & ~idx,i);
    et(oldidx & ~idx) = i;
    P(:,~idx,i) = NaN;
    oldidx = idx;
end

end

