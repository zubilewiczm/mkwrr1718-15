function [ u,t,x,y,n,m ] = CrankNicolson0( h,lambda,u0,ts,xs,ys)
%CrankNicolson0 Solves 2D Fokker-Planck Equation (that is corelated 
%with simple diffusion equation without a drift) u_t=1/2*div(Du)   
%with zero Dirichlet boundary condition, using Crank Nicolson
%finite difference method. 
%Adnotation: Because it is two dimensional equation, this
%problem is solved by using an operator splitting method called  
%ADI - Alternating direction implicit method
%   Args:
%     h         Spatial Step.
%     lambda    Constant that with given formula lambda=k/h describes 
%                 time step.
%     u0        Initial condition.
%     ts        A two-element vector representing the time interval.
%     xs        A two-element vector representing the spatial interval  
%                 on x-axis.
%     ys        A two-element vector representing the spatial interval  
%                 on y-axis.
%     boundary  String, must be either 'Dirichlet' or 'Neumann' - depends 
%                 on what boundary condition user wish to use. 
%
%   Output:
%     u         Matrix m x l x N representing solution in point (x,y) 
%                 at time t.
%     t         Time mesh.
%     x,y       Grid edges.
%     N         Length of vector t.
%     m         Length of vector x and y.

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

a = -k/(4*h.^2).*ones(m-2,1);
b = (1+k/2*h.^2).*ones(m-2,1);
c = a;
A = [0, 1, 0;
    a, b, c;
    0, 1, 0];
rhsx = zeros(1,m);
rhsy = zeros(1,l);

for n=1:N-1
    u_half(1,:) = zeros(1,l);
    for i=2:m-1
        rhsy(1) = 0;
        for j=2:l-1
            rhsy(j) = u(i,j,n) +...
                k/(4*h.^2).*(u(i,j+1,n)-2.*u(i,j,n)+u(i,j-1,n));
        end
        rhsy(l) = 0;
        u_half(i,:) = thomas(A,rhsy);
    end
    u_half(m,:) = zeros(1,l);
    
    u(:,1,n+1) = zeros(1,m);
    for j=2:l-1
        rhsx(1) = 0;
        for i=2:m-1
            rhsx(i) = u_half(i,j) ...
                + k/(4*h.^2).*(u_half(i+1,j)-2.*u_half(i,j)+u_half(i-1,j));
        end
        rhsx(m) = 0;
        u(:,j,n+1) = thomas(A,rhsx);
    end
    u(:,l,n+1) = zeros(1,m);
end