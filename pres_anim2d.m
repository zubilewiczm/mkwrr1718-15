function [ ax ] = pres_anim2d( ts, xs, tps, maxfps )
%PRES_ANIM2D Animates 2D brownian motion.
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
maxp = max(xs,[],2)*1.2;
minp = min(xs,[],2)*1.2;
plt  = plot(ax, 0,0, 0,0);
xlim(ax, [minp(1),maxp(1)]); ylim(ax, [minp(2),maxp(2)]);
daspect(ax, [1,1,1]);

plt(1).Marker = 'None';
plt(1).MarkerIndices = 1;
plt(1).MarkerEdgeColor = [1.0, 0.0, 0.0];
plt(1).MarkerSize = 12;
plt(2).Marker = 'None';
plt(2).MarkerIndices = length(plt(2).XData);
plt(2).MarkerEdgeColor = [1.0, 0.0, 0.0];
plt(2).MarkerSize = 12;

    function draw_frame(i, frames)
        if i == 1
            plt(1).Marker = 'None';
        else
            plt(1).Marker = 'o';
        end
        L = length(frames);
        fj   = frames(min(i, L));
        fjm1 = frames(max(i-1, 1));
        plt(1).XData = xs(1,1:fjm1);
        plt(1).YData = xs(2,1:fjm1);
        plt(1).Color = [0.3, 0.6, 1.0];
        plt(2).XData = xs(1,fjm1:fj);
        plt(2).YData = xs(2,fjm1:fj);
        plt(2).Color = [0.0, 0.4, 0.8];
        title(ax, sprintf('t=%.3f', ts(fj)));
        if i == L
            plt(2).Marker = 'x';
        else
            plt(2).Marker = 'None';
        end
    end

if nargin < 4
    pres_ui(@draw_frame, length(ts), ts(end)-ts(1), tps);
else
    pres_ui(@draw_frame, length(ts), ts(end)-ts(1), tps, maxfps);
end

end

