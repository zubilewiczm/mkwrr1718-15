function [ xs2, ts2 ] = brownian_refine2( xs, ts, b, r, sigma )
%BROWNIAN_REFINE2 Performs refinement of a Brownian trajectory on given set
% of time coordinates using the Brownian bridge.
%   Args:
%     xs        Positions of particle in time `ts`.
%     ts        Time mesh.
%     b         Boolean vector of the same length as `ts`. If i(j) is true,
%                 xs is refined on the interval [ts(j), ts(j+1)] with
%                 refinement ratio r(j).
%     r         Vector of refinement ratios (target step size
%                 / previous step size) of the same length as `ts`.
%                 Roughly equal to the number of subdivisions of a single
%                 step.
%     sigma     Variance matrix of position at time 1. It must be symmetric
%                 and positive-definite.
%
%   Returns:
%     xs2       Positions of particle in time `ts2`.
%     ts2       Refined time mesh.

d = size(xs,1);
if nargin < 5
    sigma = eye(d);
end

i = find(b);
k = length(i);
if isscalar(r)
    r = r*ones(1,length(k));
else
    r = r(i);
end

tt = cell(1,k);
xt = cell(1,k);

for j=1:k
    h = (ts(j+1)-ts(j))/r(j);
    [xt{j}, tt{j}] = brownian_bridge(...
        xs(:,i(j):i(j)+1), ts(i(j):i(j)+1), h, sigma);
end
sizes = cellfun(@(a) size(a,2), tt);
pts = sum(sizes)-2*k+length(ts);
ts2 = zeros(1,pts);
xs2 = zeros(d,pts);

s = 1;
p = 1; % p - i(j-1) == 1, i(0) = 0
for j = 1:k
    next = s-p+i(j);
    ts2(s:next-1)               = ts(p:i(j)-1);
    ts2(next:next+sizes(j)-2)   = tt{j}(1:end-1);
    xs2(:,s:next-1)             = xs(:,p:i(j)-1);
    xs2(:,next:next+sizes(j)-2) = xt{j}(:,1:end-1);
    s = next+sizes(j)-1;
    p = i(j)+1;
end
ts2(s:end)   = ts(p:end);
xs2(:,s:end) = xs(:,p:end);

end

