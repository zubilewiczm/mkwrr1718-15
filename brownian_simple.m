function [ xs, ts ] = brownian_simple( x0, T, h, sigma )
%BROWNIAN_SIMPLE Generate a trajectory of the Wiener process.
%   Args:
%     x0        Starting point, represented as a column vector of its 
%                 coordinates. Its length determines the dimension
%                 of the ambient space.
%     T         Duration of motion.
%     h         Time step.
%     sigma     Variance matrix of position at time 1. It must be symmetric
%                 and positive-definite.
%
%   Returns:
%     ts        Time mesh.
%     xs        Positions of particle in time `ts`.

d = length(x0);
if nargin < 4
    sigma = eye(d);
end
ts = 0:h:T;
hs = sqrt(h);
n  = length(ts);
A  = chol(sigma,'lower'); % A*A' == sigma
rv = hs*A*randn(d,n-1);
xs = cumsum([x0,rv],2);

% ...to make ts(end) == T.
if abs(ts(end)-T) >= 1e-10
    ts = [ts, T];
    h2 = sqrt(T-ts(end));
    r2 = h2*A*randn(d,1);
    xs = [xs, xs(:,end)+r2];
end

end

