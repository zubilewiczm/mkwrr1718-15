function [ P, t ] = diffuse_dirichlet( P0, bdr, T, k )
%DIFFUSE_DIRICHLET Diffuse Brownian particles inside a region in R^2
% with absorbing boundary.
%   Args:
%     P0        A dxn matrix of initial positions of particles.
%     bdr       Boolean function returning 1 if a point lies inside
%                 the region.
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
b = squeeze(bdr(P));
for i=1:n
    idx = find(~b(i,:),1);
    if ~isempty(idx)
        for j=idx+1:l
            P(:,i,j) = P(:,i,idx); % revert update
        end
    end
end
