function [ hst, hst2, x, y ] = compare_initial_cond( pdf, n, area, h )
%COMPARE_INITIAL_COND Outputs histograms over the same mesh for comparison
% of sampling methods with exact results.
%   Args:
%     pdf       Function handle representing d-dimensional probability
%                 density function of a random variable.
%     n         Number of samples.
%     area      2x2 matrix representing sampling area as a cartesian
%                 product of its rows treated as closed intervals.
%     h         Positive scalar -- mesh size.
%
%   Returns:
%     hst       Histogram of samples drawn using discrete approximation.
%     hst2      Histogram of samples drawn using the rejection method.
%     x,y       Grid edges.


% TODO: Batches for less memory footprint.
x = area(1,1):h:area(1,2);
y = area(2,1):h:area(2,2);
[X,Y] = ndgrid(x,y);
m = max(reshape(pdf(X,Y),1,[]));

samp1 = sample_discrete(pdf,n,area,h);
samp2 = sample_rejection(pdf,n,area,1.25*m);
hst   = histcounts2(samp1(1,:), samp1(2,:),...
    x, y, 'Normalization', 'pdf');
hst2  = histcounts2(samp2(1,:), samp2(2,:),...
    x, y, 'Normalization', 'pdf');

end

