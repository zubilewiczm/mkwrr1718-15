function [ P, t ] = diffuse_nobd( P0, T, k )
%DIFFUSE_NOBD Diffuse Brownian particles in R^n.
%   Args:
%     P0        A dxn matrix of initial positions of particles.
%     T         Time interval.
%     k         Time step.
%
%   Returns:
%     P         Updated positions wrt. time represented as dxnxl tensor,
%                 where l is the length of `t`.
%     t         Vector of time nodes.

t = 0:k:T;
l = length(t);
[d,n] = size(P0);
P = zeros(d,n,l);

for i=1:n
    P(:,i,:) = brownian_simple(P0(:,i),T,k);
end

end

