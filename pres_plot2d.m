function [ ax ] = pres_plot2d( ~, xs )
%PRES_PLOT2D Plots 2D brownian motion.
%   Args:
%     ts        Time mesh.
%     xs        Positions of particle in time `ts`.
%
%   Returns:
%     ax        Axes object on which the animation is drawn.

ax   = gca();
maxp = max(xs,[],2)*1.2;
minp = min(xs,[],2)*1.2;
plt  = plot(ax, xs(1,:),xs(2,:), xs(1,end),xs(2,end));
plt(1).Marker = 'o';
plt(1).MarkerIndices = 1;
plt(1).MarkerEdgeColor = [1.0, 0.0, 0.0];
plt(1).MarkerSize = 12;
plt(1).Color = [0.0, 0.4, 0.8];
plt(2).Marker = 'x';
plt(2).MarkerIndices = 1;
plt(2).MarkerEdgeColor = [1.0, 0.0, 0.0];
plt(2).MarkerSize = 12;
xlim(ax, [minp(1),maxp(1)]); ylim(ax, [minp(2),maxp(2)]);
daspect(ax, [1,1,1]);

end

