function [ xs2, ts2 ] = brownian_refine( xs, ts, i, r, sigma )
%BROWNIAN_REFINE Performs refinement of a Brownian trajectory between two
% time coordinates using the Brownian bridge.
%   Args:
%     xs        Positions of particle in time `ts`.
%     ts        Time mesh.
%     i         Index of `ts`. Refinement is done on [ts(i),ts(i+1)].
%     r         Ratio of refinement (target step size/previous step size).
%                 Roughly equals the number of subdivisions of a single
%                 step.
%     sigma     Variance matrix of position at time 1. It must be symmetric
%                 and positive-definite.
%
%   Returns:
%     xs2       Positions of particle in time `ts2`.
%     ts2       Refined ime mesh.

d = size(xs,1);
if nargin < 5
    sigma = eye(d);
end

h = (ts(i+1)-ts(i))/r;
[xt, tt] = brownian_bridge(xs(:,i:i+1), ts(:,i:i+1), h, sigma);
ts2 = [ts(1:i),tt,ts(i+1:end)];
xs2 = [xs(:,1:i),xt,xs(:,i+1:end)];

end

