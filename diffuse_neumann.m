function [ P, t ] = diffuse_neumann( P0, area, T, k )
%DIFFUSE_NEUMANN Diffuse Brownian particles inside a rectangle in R^2
% with reflecting boundary.
%   Args:
%     P0        A 2xn matrix of initial positions of particles.
%     area      A 2x2 matrix, rows of which represent bounds of the
%                 rectangle.
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

da = diff(area,1,2);
refl = @(x) area(:,2)-abs(mod(x-area(:,1), 2.*da) - da);

P = zeros(d,n,l);
for i=1:n
    P(:,i,:) = brownian_simple(P0(:,i),T,k);
    bds = squeeze(all(area(:,1) < P(:,i,:) & P(:,i,:) < area(:,2),1));
    while ~all(bds)
        s = find(~bds,1);
        v = P(:,i,s);
        w = refl(v);
        P(:,i,s:end) = P(:,i,s:end) + (w-v);
        bds = all(area(:,1) < P(:,i,:) & P(:,i,:) < area(:,2),1);
    end
end

end

