function anim_pdf( u, x, y, t, zs, fps )
%ANIM_PDF Animates the solution.
%   Args:
%     u         Solution.
%     x,y,t     Coordinates of mesh nodes.
%     zs        Z-axis limits (default: semiautomatic).
%     fps       Frames per second (default: 10).

if nargin < 5
    zs = [0,inf];
end
if nargin < 6
    fps = 10;
end

[X,Y] = ndgrid(x,y);
ax = gca();
v = 3;
for i=1:length(t)
    surf(ax, Y,X, u(:,:,i));
    view(ax,v);
    zlim(zs);
    title(ax, sprintf('Solution t=%.3f', t(i)));
    pause(1/fps);
    [a,b] = view(ax);
    v = [a,b];
end

