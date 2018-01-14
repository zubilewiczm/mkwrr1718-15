function [ ax ] = pres_anim1d( ts, xs, tps, maxfps )
%PRES_ANIM1D Animates 1D brownian motion.
%   Args:
%     ts        Time mesh.
%     xs        Positions of particle in time `ts`.
%     tps       Animation rate -- time units per second.
%     maxfps    Framerate limit (default: 30).
%
%   Returns:
%     ax        Axes object on which the animation is drawn.

fig  = gcf(); clf(fig);
ax   = gca();
maxx = max(abs(xs))*1.2;
plt  = plot(ax, 0,0, 0,0);
xlim(ax, [ts(1),ts(end)]); ylim(ax, [-maxx,maxx]);
xlabel(ax, 't'); ylabel(ax, 'B_t');
daspect(ax, [1,1,1]);

    function draw_frame(i, frames)
        L    = length(frames);
        fj   = frames(min(i, L));
        fjm1 = frames(max(i-1, 1));
        plt(1).XData = ts(1:fjm1);
        plt(1).YData = xs(1:fjm1);
        plt(1).Color = [0.3, 0.6, 1.0];
        plt(2).XData = ts(fjm1:fj);
        plt(2).YData = xs(fjm1:fj);
        plt(2).Color = [0.0, 0.4, 0.8];
    end

if nargin < 4
    pres_ui(@draw_frame, length(ts), ts(end)-ts(1), tps);
else
    pres_ui(@draw_frame, length(ts), ts(end)-ts(1), tps, maxfps);
end

end

