function [ hst, hst2, x, y ] = compare_initial_cond( pdf, n, area, h )
%COMPARE_INITIAL_COND Outputs histograms over the same mesh for comparison
% of sampling methods with exact results.
%   Args:
%     pdf       Function handle representing d-dimensional probability
%                 density function of a random variable.
%     n         Number of samples.
%     area      2x2 matrix representing sampling area as a cartesian
%                 product of its rows treated as closed intervals.
%     h         Positive scalar or 2-element vector -- mesh size.
%
%   Returns:
%     hst       Histogram of samples drawn using discrete approximation.
%     hst2      Histogram of samples drawn using the rejection method.
%     x,y       Grid edges.

if isscalar(h)
    h = [h,h];
end

x = area(1,1):h(1):area(1,2);
y = area(2,1):h(2):area(2,2);
[X,Y] = ndgrid(x,y);
batchsize = 2^16;
m = max(reshape(pdf(X,Y),1,[]));

hst  = zeros(length(x)-1, length(y)-1);
hst2 = hst;

while n>0
    k = min(batchsize, n);
    samp1 = sample_discrete(pdf,k,area,h);
    samp2 = sample_rejection(pdf,k,area,1.25*m);
    hst   = hst  + histcounts2(samp1(1,:), samp1(2,:),...
        x, y, 'Normalization', 'countdensity'); % 'countdensity'
    hst2  = hst2 + histcounts2(samp2(1,:), samp2(2,:),...
        x, y, 'Normalization', 'countdensity'); % After that normalize.
    n = n - k;
end

hst  = hst / sum(reshape(hst,1,[]));
hst2 = hst2 / sum(reshape(hst2,1,[]));

end

