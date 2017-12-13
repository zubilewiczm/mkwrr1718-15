function [ samp ] = sample_discrete( pdf, n, area, h )
%SAMPLE_DISCRETE Samples the pdf by discretizing associated r.v. and using
% its inverse cfd. Warning: memory intensive!
%   Args:
%     pdf       Function handle representing d-dimensional probability
%                 density function of a random variable.
%     n         Number of samples.
%     area      dx2 matrix representing sampling area as a cartesian
%                 product of its rows treated as closed intervals.
%     h         Positive scalar or vector of length d. Area is split into
%                 a mesh of grid size h.
%
%   Returns:
%     samp      Vector of n samples drawn from the approximate 
%                 distribution.

% setup
d = size(area,1);
g = cell(1,d);
if isscalar(h)
    h = ones(1,d)*h;
end

% discretize
s = zeros(1,d);
for i=1:d
    g{i} = area(i,1):h(i):area(i,2);
    g{i} = (g{i}(1:end-1) + g{i}(2:end))/2;
    s(i) = length(g{i});
end
[discr{1:d}] = ndgrid(g{:});

% pmf
pmf = pdf(discr{:});
if any(pmf<0)
    fails = find(pmf<0, 1);
    loc = sprintf('%.3f, ', getloc(fails, discr));
    loc = loc(1:end-2);
    error('PDF cannot be less than zero! (%s)',loc);
end

% sample
pmf = reshape(pmf,1,[]);
cdf = cumsum([0 pmf]);
cdf = cdf/cdf(end);
r   = rand(1,n);
[~, ~, ids] = histcounts(r, cdf);

% samp = zeros(d,n);
% for i=1:d
%     samp(i,:) = discr{i}(ids);
% end
samp = getloc(ids, discr);
samp = samp + (rand(d,n)-1/2).*h';

end

function locs = getloc(ids, discr)
    d = length(discr);
    n = length(ids);
    locs = zeros(d,n);
    for i=1:d
        locs(i,:) = discr{i}(ids);
    end
end
