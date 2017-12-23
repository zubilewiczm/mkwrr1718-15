function [ u,t,x,y,N,m ] = FTCSDIFF( h,lambda,u0,g,ts,xs,ys, boundary)
%FTCSDIFF Solves 2D Fokker-Planck Equation u_t(that is corelated with simple
%diffusion equation without a drift) u_t=1/2*div(Du) - using common 
%finite difference method called Forward-Time Central-Space.
%   Args:
%     h         Spatial Step.
%     lambda    Constant that with given formula lambda=k/h describes 
%                 time step.
%     u0        Initial condition.
%     g         Boundary condition.
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

[X,Y] = ndgrid(x,y);
u(:,:,1)=u0(X,Y);
    %Stability: dt<=dx^2*dy^2/(2*D*((dx^2)+(1/dy^2)))
    if k<= h^2/2                             %Remember: D=1/2, dt=k, dx=dy=h;  
    else 
    fprintf('Error, the stability condition is not met - change parameters\n')
    return
    end

for n=1:N-1
    for i=2:m-1
        for j=2:l-1
            u(i,j,n+1)=u(i,j,n)+k/2.*(...
                (u(i-1,j,n)-2*u(i,j,n)+u(i+1,j,n))./h^2 + ...
                (u(i,j-1,n)-2*u(i,j,n)+u(i,j+1,n))./h^2); 
        end
    end
    if strcmp(boundary, 'Dirichlet')

        u(1,:,n+1)=g(x(1),y,t(n+1));
        u(m,:,n+1)=g(x(m),y,t(n+1));

        u(:,1,n+1)=g(x,y(1),t(n+1));
        u(:,l,n+1)=g(x,y(l),t(n+1));

    elseif strcmp(boundary, 'Neumann')

        u(1,:,n+1)=u(2,:,n+1)+h*g(x(1),y,t(n+1));
        u(m,:,n+1)=u(m-1,:,n+1)+h*g(x(m),y,t(n+1));
    
        u(:,1,n+1)=u(:,2,n+1)+h*g(x,y(1),t(n+1)).';
        u(:,l,n+1)=u(:,l-1,n+1)+h*g(x,y(l),t(n+1)).';

    end

end

end

