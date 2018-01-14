function [ ax2 ] = ...
    pres_anim_error( ts, xs, ys, u, sol, area, titl, tps, maxfps )
%PRES_ANIM_ERROR Animate error of solution to diffusion problem
% with respect to exact solution.
%   Args:
%     ts        Time mesh.
%     xs,ys     Spatial meshes.
%     u         Approxiamte solution.
%     sol       Exact solution (function handle or another approximate
%                 solution).
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
fun = isa(sol, 'function_handle');
if size(area,1) == 2
    area = [area;-0.1,0.1];
end

fig  = gcf(); clf(fig);
ax2   = gca();
if fun
    plt2 = surf(ax2, X, Y, u(:,:,1) - sol(X,Y,ts(1)));
else
    plt2 = surf(ax2, X, Y, u(:,:,1) - sol(:,:,1));
end
xlim(ax2, area(1,:)); ylim(ax2, area(2,:)); zlim(ax2,area(3,:));

    function draw_frame(i, frames)
        if fun
            plt2.ZData = u(:,:,frames(i)) - sol(X,Y,ts(frames(i)));
        else
            plt2.ZData = u(:,:,frames(i)) - sol(:,:,frames(i));
        end
        title(ax2, sprintf('%s t=%.3f', titl, ts(frames(i))));
    end

if nargin < 9
    pres_ui(@draw_frame, length(ts), ts(end)-ts(1), tps);
else
    pres_ui(@draw_frame, length(ts), ts(end)-ts(1), tps, maxfps);
end


end

