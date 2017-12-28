% parameters
T      = 1;
area   = [0,1; -1,1];
h      = [0.1, 0.04];
lambda = 0.004;
k      = 0.004;
n      = [1e2, 1e4, 1e6];

u0 = struct('fun', @(x,y) x.*(y-0.2).^2, 'supp', area);
rg = @(x) all((area(:,1) <= x) & (x <= area(:,2)), 1);

% results
us = cell(2,3);
uf = cell(2,1);
uc = cell(2,1);
x = cell(2,1);
y = cell(2,1);

[us{1,1},x{1},y{1},t]= solve_st(u0,rg,T,area,h(1),k,n(1),'Neumann');
[us{2,1},x{2},y{2}]  = solve_st(u0,rg,T,area,h(2),k,n(1),'Neumann');
[us{1,2}]            = solve_st(u0,rg,T,area,h(1),k,n(2),'Neumann');
[us{2,2}]            = solve_st(u0,rg,T,area,h(2),k,n(2),'Neumann');
[us{1,3}]            = solve_st(u0,rg,T,area,h(1),k,n(3),'Neumann');
[us{2,3}]            = solve_st(u0,rg,T,area,h(2),k,n(3),'Neumann');

tf = cell(2,1); tc = cell(2,1);
xf = cell(2,1); xc = cell(2,1);
yf = cell(2,1); yc = cell(2,1);
[ uf{1}, tf{1}, xf{1}, yf{1} ] = ...
    solve_ftcs(h(1),lambda,u0.fun,[0,T],area(1,:),area(2,:),'Neumann');
[ uf{2}, tf{2}, xf{2}, yf{2} ] = ...
    solve_ftcs(h(2),lambda,u0.fun,[0,T],area(1,:),area(2,:),'Neumann');
[ uc{1}, tc{1}, xc{1}, yc{1} ] = ...
    solve_cn(h(1),lambda,u0.fun,[0,T],area(1,:),area(2,:),'Neumann');
[ uc{2}, tc{2}, xc{2}, yc{2} ] = ...
    solve_cn(h(2),lambda,u0.fun,[0,T],area(1,:),area(2,:),'Neumann');
% finding common times (approximate)
rat = k./(lambda.*h);
idx = round([1:rat(1):length(tf{1}); 1:rat(2):length(tf{2})]);

% Initial cond
[X,Y] = ndgrid(x{2},y{2});
figure(1);
surf(X,Y,u0.fun(X,Y));

% Errors
l = 1;
[X,Y] = ndgrid(x{l},y{l});
figure(2); ax = gca(); v = 3;
for i=1:length(t)
    surf(X,Y,us{l,3}(:,:,i)-uf{l}(:,:,idx(l,i)));
    title(sprintf('FTCS error h=%.2f t=%.2f', h(l), t(i)));
    view(ax,v);
    pause(0.04);
    [a,b] = view(ax);
    v = [a,b];
end

l = 2;
[X,Y] = ndgrid(x{l},y{l});
figure(3); ax = gca(); v = 3;
for i=1:length(t)
    surf(X,Y,us{l,3}(:,:,i)-uc{l}(:,:,idx(l,i)));
    title(sprintf('Crank-Nicolson error h=%.2f t=%.2f', h(l), t(i)));
    view(ax,v);
    pause(0.04);
    [a,b] = view(ax);
    v = [a,b];
end

% Plots
l = 1;
[X,Y] = ndgrid(x{l},y{l});
figure(4);
for i=1:length(t)
    subplot(1,2,1);
    surf(X,Y,uf{l}(:,:,idx(l,i)));
    title(sprintf('FTCS h=%.2f t=%.2f', h(l), tf{l}(idx(l,i)) ));
    
    subplot(1,2,2);
    surf(X,Y,us{l,3}(:,:,i));
    title(sprintf('Stochastic t=%.2f', t(i)));
    pause(0.04);
end

l = 2;
[X,Y] = ndgrid(x{l},y{l});
figure(5);
for i=1:length(t)
    subplot(1,2,1);
    surf(X,Y,uc{l}(:,:,idx(l,i)));
    title(sprintf('Crank-Nicolson h=%.2f t=%.2f', h(l), tc{l}(idx(l,i)) ));
    
    subplot(1,2,2);
    surf(X,Y,us{l,3}(:,:,i));
    title(sprintf('Stochastic t=%.2f', t(i)));
    pause(0.04);
end

% Max errors at
terr = [0.02, 0.4, 0.9];
usf_err = zeros(2, 3, length(terr));
usc_err = zeros(2, 3, length(terr));

for j=1:2
    for l=1:length(terr)
        for i=1:3
            usf_err(j,i,l) = ...
                max(reshape(abs(...
                    us{j,i}(:,:,t==terr(l))-uf{j}(:,:,tf{j}==terr(l))...
                ),1,[]));
            usc_err(j,i,l) = ...
                max(reshape(abs(...
                    us{j,i}(:,:,t==terr(l))-uc{j}(:,:,tc{j}==terr(l))...
                ),1,[]));
        end
    end
end

% l1 norms
us_l1 = cell(2,3);
uf_l1 = cell(2,1);
uc_l1 = cell(2,1);

for j=1:2
    uf_l1{j} = squeeze(sum(sum(uf{j}.*(h(j).^2) )));
    uc_l1{j} = squeeze(sum(sum(uc{j}.*(h(j).^2) )));
    for l=1:3
        us_l1{j,l} = squeeze(sum(sum(us{j,l}.*(h(j).^2) )));
    end
end
l=1;
figure(6);
plot(tc{l}, uc_l1{l}, t, us_l1{l,3});
legend('Crank-Nicolson', 'Stochastic n=1e+6');
title(sprintf('L1-norm h=%.2f',h(l)));

% norm drift
l=1;
figure(7);
plot(tc{l}, uc_l1{l});
legend('Crank-Nicolson');
title(sprintf('L1-norm h=%.2f',h(l)));

% maxima
usf_err2 = cell(2,3);
usc_err2 = cell(2,3);
for j=1:2
    for i=1:3
        usf_err2{j,i} = ...
            squeeze(max(max(...
            abs(us{j,i}-uf{j}(:,:,idx(j,:)))...
            )));
        usc_err2{j,i} = ...
            squeeze(max(max(...
            abs(us{j,i}-uc{j}(:,:,idx(j,:)))...
            )));
    end
end
figure(8);
subplot(1,2,1); l=1;
plot(t, usf_err2{l,3}, t, usc_err2{l,3},...
    t, movmean((usf_err2{l,3}+usc_err2{l,3})./2, 10), 'g');
title(sprintf('Max errors h=%.2f', h(l)));
legend('FTCS','Crank-Nicolson','Filtered');
subplot(1,2,2); l=2;
plot(t, usf_err2{l,3}, t, usc_err2{l,3},...
    t, movmean((usf_err2{l,3}+usc_err2{l,3})./2, 10), 'g');
title(sprintf('Max errors h=%.2f', h(l)));
legend('FTCS','Crank-Nicolson','Filtered');