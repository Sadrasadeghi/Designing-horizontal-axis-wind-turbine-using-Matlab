function [cl,cd]=interp_Radius2(cl1,cd1,cl2,cd2,r,r1,r2)
cl=cl1+(cl2-cl1)/(r2-r1)*(r-r1);
cd=cd1+(cd2-cd1)/(r2-r1)*(r-r1);
end