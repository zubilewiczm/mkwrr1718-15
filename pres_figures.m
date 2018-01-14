load('Pres.mat');

%% Brownian
figure(1);
pres_anim1d(tb1,xb1,0.5,30);

figure(2);
pres_anim2d(tb2,xb2,0.5,30);

%% Nobd delta_0
sol = @(x,y,t) 1/(2*pi*t) .* exp(-(x.^2+y.^2)/(2*t));
figure(3);
pres_anim_nobd_tr(t_nbd_tr,P_nbd,area,0.1);

figure(4);
% pres_anim_sol(t_nbd, x_nbd, y_nbd, u_nbd, sol, area, 0.1);
pres_anim_single(t_nbd, x_nbd, y_nbd, u_nbd, area, 0.1);

figure(5);
pres_anim_error(t_nbd(2:end), x_nbd, y_nbd, u_nbd, sol, area, 0.1);

%% Dirichlet
figure(6);
pres_anim_tr(td_tr, Pd, rg, area2, 0.01);

figure(7);
[X,Y] = ndgrid(xd,yd);
surf(X,Y,fun_d(X,Y));

s_ud = size(ud,3);
s_ud_cn = size(ud_cn,3);
skip_1 = (s_ud_cn-1)/(s_ud-1);

figure(8);
% pres_anim_sol(td_cn, xd, yd, ud, ud_cn(:,:,1:skip_1:end), area_d, 0.05);
pres_anim_single(td, xd, yd, ud, [area_d;0,inf], 'Ruchy Browna', 0.05);

figure(9);
pres_anim_single(td_cn, xd, yd, ud_cn, [area_d;0,inf], 'Crank-Nicolson', 0.05);

figure(10);
pres_anim_error(td, xd, yd, ud, ud_cn(:,:,1:skip_1:end), area_d, 'Błąd', 0.05);

%% Neumann
figure(11);
pres_anim_tr(tn_tr, Pn, rg, area2, 0.01);

s_un = size(un,3);
s_un_cn = size(un_cn,3);
skip_1 = (s_un_cn-1)/(s_un-1);

figure(12);
% pres_anim_sol(tn_cn, xn, yn, un, un_cn(:,:,1:skip_1:end), area_n, 0.05);
pres_anim_single(tn, xn, yn, un, area_n, 'Ruchy Browna', 0.05);

figure(13);
pres_anim_single(tn_cn, xn, yn, un_cn, area_n, 'Crank-Nicolson', 0.05);

figure(14);
pres_anim_error(tn, xn, yn, un, un_cn(:,:,1:skip_1:end), area_n, 'Błąd', 0.05);