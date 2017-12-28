% parameters
T      = 1;
area   = [-5,5; -5,5];
h      = [0.25, 0.1];
lambda = 0.004;
k      = 0.004;
n      = [1e2, 1e4, 1e6];

% results
us = cell(2,3);
uf = cell(2,1);
uc = cell(2,1);
x = cell(2,1);
y = cell(2,1);

[ us{1,1}, x{1}, y{1}, t ] = solve_st([],[],T,area,h(1),k,n(1),'None');
[ us{2,1}, x{2}, y{2} ]    = solve_st([],[],T,area,h(2),k,n(1),'None');
[ us{1,2} ]                = solve_st([],[],T,area,h(1),k,n(2),'None');
[ us{2,2} ]                = solve_st([],[],T,area,h(2),k,n(2),'None');
[ us{1,3} ]                = solve_st([],[],T,area,h(1),k,n(3),'None');
[ us{2,3} ]                = solve_st([],[],T,area,h(2),k,n(3),'None');

T0 = t(2);
u  = @(x,y,t) 1./(2*t*pi).*exp(-(x.^2+y.^2)./(2*t));
u0 = @(x,y) u(x,y,T0);
tf = cell(2,1); tc = cell(2,1);
xf = cell(2,1); xc = cell(2,1);
yf = cell(2,1); yc = cell(2,1);
[ uf{1}, tf{1}, xf{1}, yf{1} ] = ...
    solve_ftcs(h(1),lambda,u0,[T0,T],area(1,:),area(2,:),'Dirichlet');
[ uf{2}, tf{2}, xf{2}, yf{2} ] = ...
    solve_ftcs(h(2),lambda,u0,[T0,T],area(1,:),area(2,:),'Dirichlet');
[ uc{1}, tc{1}, xc{1}, yc{1} ] = ...
    solve_cn(h(1),lambda,u0,[T0,T],area(1,:),area(2,:),'Dirichlet');
[ uc{2}, tc{2}, xc{2}, yc{2} ] = ...
    solve_cn(h(2),lambda,u0,[T0,T],area(1,:),area(2,:),'Dirichlet');
% finding common times (approximate)
rat = k./(lambda.*h);
idx = round([2:1:251; 1:rat(1):length(tf{1}); 1:rat(2):length(tf{2})]);

[X,Y] = ndgrid(xf{2},yf{2});
figure(1);
surf(X,Y,u0(X,Y));

j=1;
[X,Y] = ndgrid(xf{j},yf{j});
figure(2); ax = gca(); v = 3;
for i=1:length(tf{j})
    surf(X,Y,uf{j}(:,:,i)-u(X,Y,tf{j}(i)));
    title(sprintf('FTCS error h=%.2f t=%.2f', h(j), tf{j}(i)));
    view(ax,v);
    pause(0.01);
    [a,b] = view(ax);
    v = [a,b];
end

j=2;
[X,Y] = ndgrid(xc{j},yc{j});
figure(3); ax = gca(); v = 3;
for i=1:length(tc{j})
    surf(X,Y,uc{j}(:,:,i)-u(X,Y,tc{j}(i)));
    title(sprintf('Crank-Nicolson error h=%.2f t=%.2f', h(j), tf{j}(i)));
    view(ax,v);
    pause(0.01);
    [a,b] = view(ax);
    v = [a,b];
end

j=1;
[X,Y] = ndgrid(x{j},y{j});
figure(4); ax = gca(); v = 3;
for i=1:length(t)
    surf(X,Y,us{j,3}(:,:,i)-u(X,Y,t(i)));
    title(sprintf('Stochastic error h=%.2f t=%.2f', h(j), tf{j}(i)));
    view(ax,v);
    pause(0.01);
    [a,b] = view(ax);
    v = [a,b];
end

j=1;
figure(5);
[X,Y] = ndgrid(xc{j},yc{j});
[X2,Y2] = ndgrid(x{j},y{j});
terr = 0.94;
subplot(1,2,1);
surf(X,Y,uc{j}(:,:,tc{j}==terr)-u(X,Y,terr));
title(sprintf('Crank-Nicloson h=%.2f t=%.2f',h(1),terr));
subplot(1,2,2);
surf(X,Y,us{j,3}(:,:,t==terr)-u(X,Y,terr));
title(sprintf('Stochastic N=1e+06 h=%.2f t=%.2f',h(1),terr));

% Max errors at
terr = [0.02, 0.4, 0.9];
us_err = zeros(2, 3, length(terr));
uf_err = zeros(2, length(terr));
uc_err = zeros(2, length(terr));
uf_err_bd = zeros(2, length(terr));
uc_err_bd = zeros(2, length(terr));

for j=1:2
    [X,Y] = ndgrid(x{j},y{j});
    for l=1:length(terr)
        for i=1:3
            us_err(j,i,l) = ...
                max(reshape(abs(us{j,i}(:,:,t==terr(l))-u(X,Y,terr(l))),1,[]));
        end
        uf_err(j,l) = ...
            max( reshape(abs(uf{j}(:,:,tf{j}==terr(l))-u(X,Y,terr(l))),1,[]));
        uc_err(j,l) = ...
            max( reshape(abs(uc{j}(:,:,tc{j}==terr(l))-u(X,Y,terr(l))),1,[]));
        uf_err_bd(j,l) = ...
            max( abs(uf{j}(1,:,tf{j}==terr(l))-u(x{j}(1),y{j},terr(l))) );
        uc_err_bd(j,l) = ...
            max( abs(uc{j}(1,:,tc{j}==terr(l))-u(x{j}(1),y{j},terr(l))) );
    end
end

% Max errors total
us_err2 = cell(2,3);
uf_err2 = cell(2,1);
uc_err2 = cell(2,1);
terr2   = t(idx(1,:));
for j=1:2
    [X,Y] = ndgrid(x{j},y{j});
    utmp = zeros(length(x{j}), length(y{j}), length(terr2));
    for i=1:length(terr2)
        utmp(:,:,i) = u(X,Y,terr2(i));
    end
    for i=1:3
        us_err2{j,i} = ...
            squeeze(max(max(...
            abs(us{j,i}(:,:,idx(1,:))-utmp)...
            )));
    end
    uf_err2{j} = ...
        squeeze(max(max(...
            abs(uf{j}(:,:,idx(j+1,:))-utmp)...
        )));
    uc_err2{j} = ...
        squeeze(max(max(...
            abs(uc{j}(:,:,idx(j+1,:))-utmp)...
        )));
end
j=1;
figure(6);
plot(terr2, uf_err2{j}, terr2, uc_err2{j}, terr2, us_err2{j,3});
title(sprintf('Max errors h=%.2f', h(j)));
legend('FTCS','Crank-Nicolson','Stochastic');
xlabel('t'); ylabel('max u');