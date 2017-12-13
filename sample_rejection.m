function [ samp ] = sample_rejection( pdf, n, area, top )
%SAMPLE_REJECTION Samples the pdf using the rejection sampling method.
%   Args:
%     pdf       Function handle representing d-dimensional probability
%                 density function of a random variable.
%     n         Number of samples.
%     area      dx2 matrix representing sampling area as a cartesian
%                 product of its rows treated as closed intervals.
%     top       Cutoff value. After choosing a point in the sampling
%                 area take another uniformly random real from the interval
%                 [0, top]. If it's less than the pdf at the point, accept
%                 the sample; reject otherwise.
%
%   Returns:
%     samp      Vector of n samples drawn from given distribution.

draw = 2.^16;

d     = size(area,1);
delta = diff(area,1,2);

samp = zeros(d, n);
l    = 0;
while l < n
    r = rand(d+1,draw) .* [delta; top] + [area(:,1); 0];
    c = cell(1,d);
    for i=1:d
        c{i} = r(i,:);
    end
    r = r(:, pdf(c{:}) > r(d+1,:));
    s = min(size(r,2),n-l);
    samp(:,l+1:l+s) = r(1:d,1:s);
    l = l+s;
end

end

