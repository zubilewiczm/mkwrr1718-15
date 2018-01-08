function [ ax ] = pres_plot1d( ts, xs )
%PRES_PLOT1D Plots 1D brownian motion.
%   Args:
%     ts        Time mesh.
%     xs        Positions of particle in time `ts`.
%
%   Returns:
%     ax        Axes object on which the animation is drawn.

ax   = gca();
maxx = max(abs(xs))*1.2;
plot(ax, ts, xs, 'Color', [0.0, 0.4, 0.8]);
xlim(ax, [0,10]); ylim(ax, [-maxx,maxx]);
xlabel(ax, 't'); ylabel(ax, 'B_t');
daspect(ax, [1,1,1]);

end

