function [ ax ] = pres_anim_tr( ts, P, bdr, area, tps, maxfps )
%PRES_ANIM_BD_TR Animate trajectories of multiple brownian particles.
%Highlight the domain's boundary.
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

% Scatter plot
cmap = [0,   0.8, 0;
        0.2, 0,   0.8;
        0.8, 0,   0.2;
        0.2, 0.5, 0.2];
c    = interp1(cmap, linspace(1,size(cmap,1),size(P,2)));
plt  = scatter(ax, P(1,:,1) ,P(1,:,1), 12, c);
xlim(ax, area(1,:)); ylim(ax, area(2,:));
daspect(ax, [1,1,1]);

% Boundary
xs = linspace(area(1,1), area(1,2), 300);
ys = linspace(area(2,1), area(2,2), 300);
[X,Y] = ndgrid(xs, ys);
pts = [reshape(X,1,[]); reshape(Y,1,[])];
idx = bdr(pts);
pts = pts(:,idx);
idx_bd  = boundary(pts.', 1);
pts = pts(:,idx_bd);
line(pts(1,:), pts(2,:), 'Color', 'k', 'LineStyle', '--');

    function draw_frame(i, frames)
        fj = frames(i);
        inbd = bdr(P(:,:,fj));
        plt.XData = P(1,:,fj);
        plt.YData = P(2,:,fj);
        plt.CData = (c .* inbd.') + [1.0, 0.0, 0.0] .* (1-inbd.');
        title(ax, sprintf('t=%.3f', ts(fj)));
    end

if nargin < 6
    pres_ui(@draw_frame, length(ts), ts(end)-ts(1), tps);
else
    pres_ui(@draw_frame, length(ts), ts(end)-ts(1), tps, maxfps);
end


end

