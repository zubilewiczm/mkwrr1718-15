pdf  = @(x,y) x.*y.^2;
area = [0,1;-1,1];
k    = 10;
h    = (diff(area,1,2)./k).';

n = 10.^(2:9);
[itg, x, y] = compare_initial_itg(pdf,area,h);

hst  = cell(1,length(n));
hst2 = cell(1,length(n));

% Compute
for i=1:length(n)
    [hst{i}, hst2{i}] = compare_initial_cond(pdf,n(i),area,h);
end

% Display
figure(1);
x = linspace(area(1,1), area(1,2), 20);
y = linspace(area(2,1), area(2,2), 20);
[X,Y] = ndgrid(x,y);
Z = pdf(X,Y);
surf(X,Y,Z);
xlabel('x'); ylabel('y');

for i=1:length(n)
    k = n(i);
    figure(i+1);
    subplot(1,2,1);
    histogram2('XBinEdges',x,'YBinEdges',y,'BinCounts',abs(hst{i}-itg));
    title(sprintf('Discretization n=%.0e', k));
    subplot(1,2,2);
    histogram2('XBinEdges',x,'YBinEdges',y,'BinCounts',abs(hst2{i}-itg));
    title(sprintf('Rejection n=%.0e', k));
end

% Get max errors
maxerr = zeros(2,length(n));
for i=1:length(n)
    maxerr(1,i) = max(reshape(abs(hst{i}-itg),1,[]));
    maxerr(2,i) = max(reshape(abs(hst2{i}-itg),1,[]));
end