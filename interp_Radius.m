function [cl,cd,betta_inf]=interp_Radius(cl1,cd1,cl2,cd2,alpha1,alpha2,r,r1,r2)
% This function will find the values of cl and cd at different radius by
% using of the linear interpolation.wehave the cl and cd and certaint
% points so now we must calculate them for other radius(r/R)

[eff1,p1]=max(cl1./cd1);
[eff2,p2]=max(cl2./cd2);

cl=cl1(p1)+(cl2(p2)-cl1(p1))/(r2-r1)*(r-r1);
cd=cd1(p1)+(cd2(p2)-cd1(p1))/(r2-r1)*(r-r1);
betta_inf=alpha1(p1)+(alpha2(p2)-alpha1(p1))/(r2-r1)*(r-r1);

end
