function [x_t] = new_x(x,u,n)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
x_t = zeros(n, 2);
j=1;
while j<=n
    x_t(j, 1) = x(j,1) + u(2*j-1,1);
    x_t(j, 2) = x(j,2) + u(2*j, 1);
    j=j+1;
end

