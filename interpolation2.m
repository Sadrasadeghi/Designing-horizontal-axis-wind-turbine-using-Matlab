function [a]=interpolation2(a)
b=a(1,1):0.01:a(end,1);
%% The function presented here interpolates cl&cd values at different angles.
cl=spline(a(:,1),a(:,2),b);
cd=spline(a(:,1),a(:,3),b);

a=[b;cl;cd]';

end