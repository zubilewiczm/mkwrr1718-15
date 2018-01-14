function [ ax1 ] = ...
    pres_anim_single( ts, xs, ys, u, area, titl, tps, maxfps )
%PRES_ANIM_SINLGE Animate solution to diffusion problem.
%   Args:
%     ts        Time mesh.
%     xs,ys     Spatial meshes.
%     u         Approxiamte solution.
%     area      A 2x2 (or 2x3) matrix representing solution's bounding box
%                 as a cartesian product of its rows treated as closed
%                 intervals.
%     titl      Plot title.
%     tps       Animation rate -- time units per second.
%     maxfps    Framerate limit (default: 30).
%
%   Returns:
%     ax        Axes object on which the animation is drawn.

[X, Y] = ndgrid(xs,ys);
if size(area,1) == 2
    area = [area;0,inf];
end

fig  = gcf(); clf(fig);
ax1   = gca();
plt1  = surf(ax1, X, Y, u(:,:,1));
xlim(ax1, area(1,:)); ylim(ax1, area(2,:)); zlim(ax1,area(3,:));
% daspect(ax1, [1,1,1]); daspect(ax2, [1,1,1]);

    function draw_frame(i, frames)
        plt1.ZData = u(:,:,frames(i));
        title(ax1, sprintf('%s t=%.3f', titl, ts(frames(i))));
    end

if nargin < 9
    pres_ui(@draw_frame, length(ts), ts(end)-ts(1), tps);
else
    pres_ui(@draw_frame, length(ts), ts(end)-ts(1), tps, maxfps);
end


end

