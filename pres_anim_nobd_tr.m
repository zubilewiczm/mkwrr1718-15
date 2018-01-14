function [ ax ] = pres_anim_nobd_tr( ts, P, area, tps, maxfps )
%PRES_ANIM_NOBD_TR Animate trajectories of multiple brownian particles.
%   Args:
%     ts        Time mesh.
%     P         Positions of particles in time `ts`, in form of a 2xnxm
%                 tensor, where `n` is the number of particles
%                 and `m == length(ts)`.
%     area      A 2x2 matrix representing trajectories' bounding box
%                 as a cartesian product of its rows treated as closed
%                 intervals.
%     tps       Animation rate -- time units per second.
%     maxfps    Framerate limit (default: 30).
%
%   Returns:
%     ax        Axes object on which the animation is drawn.

fig  = gcf(); clf(fig);
ax   = gca();
cmap = [0,   0.8, 0;
        0.2, 0,   0.8;
        0.8, 0,   0.2;
        0.2, 0.5, 0.2];
c    = linspace(1,size(cmap,1),size(P,2));
plt  = scatter(ax, P(1,:,1) ,P(1,:,1), 12, interp1(cmap,c));
xlim(ax, area(1,:)); ylim(ax, area(2,:));
daspect(ax, [1,1,1]);

    function draw_frame(i, frames)
        fj = frames(i);
        plt.XData = P(1,:,fj);
        plt.YData = P(2,:,fj);
        title(ax, sprintf('t=%.3f', ts(fj)));
    end

if nargin < 5
    pres_ui(@draw_frame, length(ts), ts(end)-ts(1), tps);
else
    pres_ui(@draw_frame, length(ts), ts(end)-ts(1), tps, maxfps);
end


end

