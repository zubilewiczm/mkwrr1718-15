pdf  = @(x,y) x.*y.^2./5;
area = [0,1;-1,1];
h    = [0.1, 0.2];

hst  = cell(1,5);
hst2 = cell(1,5);

n = 100;
[itg, x, y] = compare_initial_itg(pdf,area,h);

for i=1:3
    k = n*100^i;
    [hst{i}, hst2{i}] = compare_initial_cond(pdf,k,area,h);
    figure(i);
    subplot(1,2,1);
    histogram2('XBinEdges',x,'YBinEdges',y,'BinCounts',abs(hst{i}-itg));
    title(sprintf('Discretization n=%d', k));
    subplot(1,2,2);
    histogram2('XBinEdges',x,'YBinEdges',y,'BinCounts',abs(hst2{i}-itg));
    title(sprintf('Rejection n=%d', k));
end