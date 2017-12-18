function [ u, x, y, t ] = solve_st_nobd( u0, T, area, h, k, n  )
%SOLVE_ST_NOBD Solves the diffusion equation on R^2 using Brownian motion.
%   Args:
%     u0        struct
%                   .fun    Function returning initial values.
%                   .supp   Support of the function.
%                   .m      (opional) Maximal value of the function.
%               or a function handle corresponding to u0.fun above.
%     T         Length of time interval.
%     area      2x2 matrix representing sampling area as a cartesian
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

hasm = 0;
if isstruct(u0)
    if isfield(u0,'m')
        m = u0.m;
        hasm = 1;
    else
        supp = u0.supp;
        u0 = u0.fun;
    end
else
    supp = area;
end

if ~hasm
    res = 1000;
    x = linspace(supp(1,1), supp(1,2), res);
    y = linspace(supp(2,1), supp(2,2), res);
    [X,Y] = ndgrid(x,y);
    m = max(reshape(u0(X,Y),1,[]));
end

itg = integral2(u0, supp(1,1), supp(1,2), supp(2,1), supp(2,2));

P0 = sample_rejection(u0,n,supp,1.25*m);
[P, t] = diffuse_nobd(P0,T,k);

[u_tmp, x, y] = normalize(P(:,:,1), itg, n, area, h);
lt = length(t); lx = length(x); ly = length(y);
u = zeros(lx, ly, lt);
u(:,:,1) = u_tmp;
for i=2:lt
    u(:,:,i) = normalize(P(:,:,i), itg, n, area, h);
end

end
