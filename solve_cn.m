function [ u,t,x,y,N,m ] = solve_cn( h,lambda,u0,ts,xs,ys,boundary )
%SOLVE_CN Solves 2D Diffusion equation withouth a drift: u_t=1/2*div(Du)   
%with zero Dirichlet/Neumann boundary conditions, using Crank Nicolson
%finite difference method. 
%
%Adnotation1: Because it is two dimensional equation, this
%problem is solved by using an operator splitting method called  
%ADI - Alternating direction implicit method
%
%Adnotation2: This function is an improvement over both functions
%CrankNicolson0N and CrankNicolson0N. There are added additional spectral
%points lying beyond our initial mesh, that helps finding second order
%approximations of boundary conditions.
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
x=xs(1)+h/2:h:xs(2)-h/2;
y=ys(1)+h/2:h:ys(2)-h/2;
t=ts(1):k:ts(2);

N=length(t);
m=length(x);
l=length(y);

u=zeros(m,l,N);
u_tmp  = zeros(m+2,l+2);
u_half = zeros(m+2,l+2);

[X,Y] = ndgrid(x,y);
u(:,:,1)=u0(X,Y);

a = -k/(2*h.^2);      ax = a.*ones(m,1); ay = a.*ones(l,1);
b = 2*(1+k/(2*h.^2)); bx = b.*ones(m,1); by = b.*ones(l,1);
c = a;                cx = c.*ones(m,1); cy = c.*ones(l,1);

if strcmp(boundary, 'Dirichlet')
    Ax = [0,  1,  1;
          ax, bx, cx;
           1, 1,  0];
    Ay = [0,  1,  1;
          ay, by, cy;
           1, 1,  0];
       
elseif strcmp(boundary, 'Neumann')
    Ax = [0,  1, -1;
          ax, bx, cx;
          -1, 1,  0];
    Ay = [0,  1, -1;
          ay, by, cy;
          -1, 1,  0];
end
rhsx = zeros(1,m+2);
rhsy = zeros(1,l+2);

for n=1:N-1
    u_tmp(2:end-1,2:end-1) = u(:,:,n);
    % approximations of order 2
    if strcmp(boundary, 'Dirichlet')
        u_tmp(1,:)   = -u_tmp(2,:);
        u_tmp(m+2,:) = -u_tmp(m+1,:);
        u_tmp(:,1)   = -u_tmp(:,2);
        u_tmp(:,l+2) = -u_tmp(:,l+1);
    elseif strcmp(boundary, 'Neumann')
        u_tmp(1,:)   = u_tmp(2,:);
        u_tmp(m+2,:) = u_tmp(m+1,:);
        u_tmp(:,1)   = u_tmp(:,2);
        u_tmp(:,l+2) = u_tmp(:,l+1);
    end
    for i=1:m
        rhsy(1) = 0;
        for j=1:l
            rhsy(j+1) = 2*u_tmp(i+1,j+1) + k/(2*h.^2).* ...
                (u_tmp(i+2,j+1)-2.*u_tmp(i+1,j+1)+u_tmp(i,j+1));
        end
        rhsy(l+2) = 0;
        u_half(i+1,:) = thomas(Ay,rhsy);
    end
    if strcmp(boundary, 'Dirichlet')
        u_half(1,:)   = -u_half(2,:);
        u_half(m+2,:) = -u_half(m+1,:);
    elseif strcmp(boundary, 'Neumann')
        u_half(1,:)   = u_half(2,:);
        u_half(m+2,:) = u_half(m+1,:);
    end
    
    for j=1:l
        rhsx(1) = 0;
        for i=1:m
            rhsx(i+1) = 2*u_half(i+1,j+1) + k/(2*h.^2).* ...
                (u_half(i+1,j+2)-2.*u_half(i+1,j+1)+u_half(i+1,j));
        end
        rhsx(m+2) = 0;
        u_tmp(:,j+1) = thomas(Ax,rhsx);
    end
    u(:,:,n+1) = u_tmp(2:end-1,2:end-1);
end
end
