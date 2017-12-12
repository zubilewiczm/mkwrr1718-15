%% 2D Brownian motion
x1 = brownian_simple([0;0],1,0.001);
figure(1);
plot(x1(1,:),x1(2,:));
daspect([1,1,1]);

%% 1D Brownian motion vs time
[x2,t2] = brownian_simple(0,1,0.03);
figure(2);
plot(t2,x2);

%% 1D Brownian bridge
[x3,t3] = brownian_bridge([-1,1],[0,1],0.001);
figure(3);
plot(t3,x3);

%% 2D Brownian bridge
x4 = brownian_bridge([-1,1;1,-1],[0,1],0.0001);
figure(4);
plot(x4(1,:),x4(2,:));
daspect([1,1,1]);

%% 1D Refinement
[x5,t5] = brownian_refine(x2,t2,15,100);
figure(5);
plot(t5,x5,t2,x2);

[x6,t6] = brownian_refine2(x2,t2,...
    (t2 > 0.3) & (t2 < 0.7),... & mod(1:length(t2),2)==0,...
    max(10,(x2+0.3)*400));
figure(6);
plot(t6,x6,t2,x2);

%% Distribution
d1 = diff(x1,1,2);
figure(7);
subplot(1,2,1);
scatter(d1(1,:), d1(2,:));
xlim([-0.15,0.15]);
ylim([-0.15,0.15]);
subplot(1,2,2);
histogram2(d1(1,:), d1(2,:), 'Normalization', 'pdf');

hold on;
I = linspace(-0.3,0.3, 30);
[X,Y] = meshgrid(I,I);
Z = 1/(2*pi*0.001)*exp(-(X.^2+Y.^2)/0.002);
surf(X,Y,Z);
alpha(0.7);
hold off;