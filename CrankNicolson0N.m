function [ u,t,x,y,n,m ] = CrankNicolson0N( h,lambda,u0,ts,xs,ys)

k=lambda*h;
x=xs(1):h:xs(2);
y=ys(1):h:ys(2);
t=ts(1):k:ts(2);

N=length(t);
m=length(x);
l=length(y);

u=zeros(m,l,N);
u_half = zeros(m,l);

[X,Y] = ndgrid(x,y);
u(:,:,1)=u0(X,Y);

a = -k/(2*h.^2);      ax = a.*ones(m-2,1); ay = a.*ones(l-2,1);
b = 2*(1+k/(2*h.^2)); bx = b.*ones(m-2,1); by = b.*ones(l-2,1);
c = a;                cx = c.*ones(m-2,1); cy = c.*ones(l-2,1);
Ax = [0, 1, -1;
      ax, bx, cx;
      -1, 1, 0];
Ay = [0, 1, -1;
      ay, by, cy;
   -  1, 1, 0];
rhsx = zeros(1,m);
rhsy = zeros(1,l);

for n=1:N-1
   for i=2:m-1
       rhsy(1) = 0;
       for j=2:l-1
           rhsy(j) = 2*u(i,j,n) +...
               k/(2*h.^2).*(u(i+1,j,n)-2.*u(i,j,n)+u(i-1,j,n));
       end
       rhsy(l) = 0;
       u_half(i,:) = thomas(Ay,rhsy);
   end
   u_half(1,:) = u_half(2,:);
   u_half(m,:) = u_half(m-1,:);
   
   for j=2:l-1
       rhsx(1) = 0;
       for i=2:m-1
           rhsx(i) = 2*u_half(i,j) ...
               + k/(2*h.^2).*(u_half(i,j+1)-2.*u_half(i,j)+u_half(i,j-1));
       end
       rhsx(m) = 0;
       u(:,j,n+1) = thomas(Ax,rhsx);
   end
   u(:,1,n+1) = u(:,2,n+1);
   u(:,l,n+1) = u(:,l-1,n+1);
%    u(1,1,n+1) = (u(2,1,n+1) + u(1,2,n+1))./2;
%    u(m,1,n+1) = (u(m-1,1,n+1) + u(m,2,n+1))./2;
%    u(1,l,n+1) = (u(2,l,n+1) + u(1,l-1,n+1))./2;
%    u(m,l,n+1) = (u(m-1,l,n+1) + u(m,l-1,n+1))./2;
end