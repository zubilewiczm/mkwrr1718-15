function [ x ] = thomas( A, y )
%THOMAS Thomas Algorithm is a simplified form of Gaussian elimination 
%that solves tridiagonal systems of equations
%   A -> Matrix nx3
%   y -> vector of n elements

    n = length(y);
    p = zeros(n-1,1);
    q = p;
    x = zeros(n,1);
    a = A(:,1);
    b = A(:,2);
    c = A(:,3);
    
    p(1) = -c(1)/b(1);
    q(1) = y(1);
    for i=2:n-1
        p(i) = -c(i)./(a(i).*p(i-1)+b(i));
        q(i) = (y(i) - a(i).*q(i-1))/(a(i).*p(i-1) + b(i));
    end
    x(n) = (y(n)-a(n).*q(n-1))/(b(n)+a(n).*p(n-1));
    for i=n-1:-1:1
        x(i) = p(i).*x(i+1)+q(i);
    end
    
end

