function [ xs, ts ] = brownian_bridge( xi, ti, h, sigma )
%BROWNIAN_BRIDGE Generate a trajectory of the Wiener process convergent to
% a given position at chosen time.
%   Args:
%     xi        Initial and final position of the particle, represented
%                 as a dx2 matrix, where d is the dimension
%                 of the ambient space.
%     ti        A two-element vector representing the time interval, during
%                 which the motion occurs.
%     h         Time step.
%     sigma     Variance matrix of position at time 1. It must be symmetric
%                 and positive-definite.
%
%   Returns:
%     xs        Positions of particle in time `ts`.
%     ts        Time mesh.

d = size(xi,1);
if nargin < 4
    sigma = eye(d);
end

[xs, ts] = brownian_simple(xi(:,1), ti(2)-ti(1), h, sigma);
xs = xs + ts/(ti(2)-ti(1)).*(xi(:,2)-xs(:,end));
ts = ts + ti(1);

end

