function [ u,t,x,y,N,m ] = solve_ftcs( h,lambda,u0,ts,xs,ys, boundary)
%SOLVE_FTCS Solves 2D diffusion equation  u_t=1/2*div(Du) with zero
%Dirichlet/Neumann boundary conditions using common 
%finite difference method called Forward-Time Central-Space.
%Function itself is any improvement over FTCSDIFF.m. Spatial grid was
%extended over h/2 step in both directions.This helps with providing  
%better numerical approximations of boundary conditions.
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

[X,Y] = ndgrid(x,y);
u(:,:,1)=u0(X,Y);
    %Stability: dt<=dx^2*dy^2/(2*D*((dx^2)+(1/dy^2)))
    if k<= h^2/2                             %Remember: D=1/2, dt=k, dx=dy=h;  
    else 
    fprintf('Error, the stability condition is not met - change parameters\n')
    return
    end

u_tmp = zeros(m+2,l+2);
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
        for j=1:l
            u(i,j,n+1)=u_tmp(i+1,j+1)+k/2.*(...
                u_tmp(i,j+1)-2*u_tmp(i+1,j+1)+u_tmp(i+2,j+1) + ...
                u_tmp(i+1,j)-2*u_tmp(i+1,j+1)+u_tmp(i+1,j+2))./h^2;
        end
    end
end

end
