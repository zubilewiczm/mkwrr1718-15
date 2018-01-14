function [ ax1, ax2 ] = ...
    pres_anim_sol( ts, xs, ys, u, sol, area, tps, maxfps )
%PRES_ANIM_SOL Animate solution to diffusion problem.
% Compare approximate solution with analytic one.
%   Args:
%     ts        Time mesh.
%     xs,ys     Spatial meshes.
%     u         Approxiamte solution.
%     sol       Exact solution (function handle or another approximate
%                 solution).
%     area      A 2x2 (or 2x3) matrix representing solution's bounding box
%                 as a cartesian product of its rows treated as closed
%                 intervals.
%     tps       Animation rate -- time units per second.
%     maxfps    Framerate limit (default: 30).
%
%   Returns:
%     ax        Axes object on which the animation is drawn.

[X, Y] = ndgrid(xs,ys);
if size(area,1) == 2
    area = [area;0,inf];
end
fun = isa(sol, 'function_handle');

fig  = gcf(); clf(fig);
ax1   = subplot(1,2,1);
ax2   = subplot(1,2,2);
plt1  = surf(ax1, X, Y, u(:,:,1));
if fun
    plt2 = surf(ax2, X, Y, u(:,:,1) - sol(X,Y,ts(1)));
else
    plt2 = surf(ax2, X, Y, u(:,:,1) - sol(:,:,1));
end
xlim(ax1, area(1,:)); ylim(ax1, area(2,:)); zlim(ax1,area(3,:));
xlim(ax2, area(1,:)); ylim(ax2, area(2,:)); zlim(ax2,[-0.25,0.25]);
% daspect(ax1, [1,1,1]); daspect(ax2, [1,1,1]);

    function draw_frame(i, frames)
        plt1.ZData = u(:,:,frames(i));
        if fun
            plt2.ZData = u(:,:,frames(i)) - sol(X,Y,ts(frames(i)));
        else
            plt2.ZData = u(:,:,frames(i)) - sol(:,:,frames(i));
        end
        title(ax1, sprintf('Ruchy Browna t=%.3f', ts(frames(i))));
        title(ax2, sprintf('Błąd t=%.3f',  ts(frames(i))));
    end

if nargin < 8
    pres_ui(@draw_frame, length(ts), ts(end)-ts(1), tps);
else
    pres_ui(@draw_frame, length(ts), ts(end)-ts(1), tps, maxfps);
end


end

