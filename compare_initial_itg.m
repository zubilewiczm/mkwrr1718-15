function [ itg, x, y ] = compare_initial_itg( pdf, area, h )
%COMPARE_INITIAL_ITG Outputs integrals over cells of a mesh, normalized
% over area.
%   Args:
%     pdf       Function handle representing d-dimensional probability
%                 density function of a random variable.
%     area      2x2 matrix representing sampling area as a cartesian
%                 product of its rows treated as closed intervals.
%     h         Positive scalar -- mesh size.
%
%   Returns:
%     itg       Integrals over cells of a mesh.
%     x,y       Grid edges.

x = area(1,1):h:area(1,2);
y = area(2,1):h:area(2,2);
[X,Y] = ndgrid(x,y);
itg = arrayfun(@(x1,x2,y1,y2) integral2(pdf, x1,x2, y1,y2),...
    X(1:end-1,1:end-1), X(2:end,1:end-1),...
    Y(1:end-1,1:end-1), Y(1:end-1,2:end));
itg = itg / sum(reshape(itg,1,[]));

end

