function [ epdf, x, y ] = normalize( P, itg, n, area, h )
%NORMALIZE Extract pdf from 2D empirical distribution.
%   Args:
%     P         A posteriori empirical distribution in the form of 
%                 a 2xN matrix of column samples.
%     itg       Normalization coefficent.
%     n         Number of inital samples.
%     area      2x2 matrix representing sampling area as a cartesian
%                 product of its rows treated as closed intervals.
%     h         Positive scalar or 2-element vector -- mesh size.
%
%   Returns:
%     epdf      Computed pdf, normalized so that the integral of epdf over
%                 the area equals the integral of pdf0 * N/n.
%     x,y       Grid edges.

if isscalar(h)
    h = [h,h];
end

x = area(1,1):h(1):area(1,2);
y = area(2,1):h(2):area(2,2);

epdf = itg .* histcounts2(P(1,:), P(2,:), x, y,...
    'Normalization', 'countdensity') ./ n;
x = movmean(x,2, 'Endpoints', 'discard');
y = movmean(y,2, 'Endpoints', 'discard');

end

