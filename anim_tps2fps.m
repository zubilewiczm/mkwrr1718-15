function [ fps, ind ] = anim_tps2fps( tps, T, N, maxfps )
%ANIM_TPS2FPS Picks indices `ind` of `1:length(ts)` in such a way, that
%`ts(ind(k:k+fps))` contains times from `ts(ind(k))` to `ts(ind(k)) + tps`.
%   This lets the user make an animation that shows `tps` time units of
%   a simulation in 1 real second.
%
%   Args:
%     tps       Rate of animation.
%     ts        Animation time points.
%     maxfps    Framerate limit (default: 30).
%
%   Returns:
%     fps       Framerate. May be less than `maxfps` if there aren't many
%                 time points in `ts`.
%     ind       Subsequence of `ts`.

if nargin < 4
    maxfps = 30;
end

expfps = tps*(N-1)/T;
fps = min(expfps, maxfps);
ind = [round(1:expfps/fps:N), N];

end

