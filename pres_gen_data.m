%% Brownian
T = 10;
h = 0.0001;
[xb1, tb1] = brownian_simple(0, T, h);
[xb2, tb2] = brownian_simple([0;0], T, h);

%% Nobd
area = [-5,5;-5,5];
T = 3;
h = 0.1;
k = 0.004;
[u_nbd,x_nbd,y_nbd,t_nbd] = ...
    solve_st([], [], T, area, h, k, 1000000, 'None');
T = 3;
k = 0.0001;
[~,~,~,t_nbd_tr,P_nbd] = ...
    solve_st([], [], T, area, h, k, 2000, 'None');

%% Dirichlet
area2 = [-1.2,1.2;-1.2,1.2];
bdr = [-1,1;-1,1];
T = 1;
h = 0.1;
k = 0.0001;
rg = @(x) all(abs(x) < 1);
[~,~,~,td_tr,Pd] = solve_st([], rg, T, bdr, h, k, 500, 'Dirichlet');

T = 1;
area_d = [0,1; -1,1];
lambda = 0.004;
k      = 0.004;
h      = 0.1;
fun_d = @(x,y) x.*(y-0.2).^2;
u0d = struct('fun', fun_d, 'supp', area_d);
rgd = @(x) all((area_d(:,1) <= x) & (x <= area_d(:,2)), 1);
[ud, xd, yd, td] = ...
    solve_st(u0d, rgd, T, area_d, h, k, 1000000, 'Dirichlet');

[ud_cn, td_cn] = solve_cn(h, lambda, fun_d, ...
    [0,T], area_d(1,:), area_d(2,:), 'Dirichlet');

%% Neumann
T = 1;
h = 0.1;
k = 0.0001;
[~,~,~,tn_tr,Pn] = solve_st([], rg, T, bdr, h, k, 1000, 'Neumann');

T = 1;
area_n = [0,1; -1,1];
lambda = 0.004;
k      = 0.004;
h      = 0.1;
fun_n = @(x,y) x.*(y-0.2).^2;
u0d = struct('fun', fun_n, 'supp', area_n);
rgd = @(x) all((area_n(:,1) <= x) & (x <= area_n(:,2)), 1);
[un, xn, yn, tn] = ...
    solve_st(u0d, rgd, T, area_n, h, k, 1000000, 'Neumann');

[un_cn, tn_cn] = solve_cn(h, lambda, fun_n, ...
    [0,T], area_n(1,:), area_n(2,:), 'Neumann');