function [CL,CD,inc_angle]=interp_Re(a10,a20,a30,a40,Re)
% Reynolds Interpolation
% cl and cd at the same incident angle will be changed when we change
% Reynolds number so here we define a function through which we can compute
% cl and cd at different Re.

inc_angle=a10(:,1); % we have constant alpha

range=1e6:1e6:4e6; 


cl=[a10(:,2),a20(:,2),a30(:,2),a40(:,2)]'; 
cd=[a10(:,3),a20(:,3),a30(:,3),a40(:,3)]'; 

if Re<range(1)
      CL=cl(1,:);
      CD=cd(1,:);
elseif Re>=range(1) && Re<range(2)
      CL=cl(1,:)+(cl(2,:)-cl(1,:))/(range(2)-range(1))*(Re-range(1));
      CD=cd(1,:)+(cd(2,:)-cd(1,:))/(range(2)-range(1))*(Re-range(1));
elseif Re>=range(2) && Re<range(3)
      CL=cl(2,:)+(cl(3,:)-cl(2,:))/(range(3)-range(2))*(Re-range(2));
      CD=cd(2,:)+(cd(3,:)-cd(2,:))/(range(3)-range(2))*(Re-range(2));
elseif Re>=range(3) && Re<range(4)
      CL=cl(3,:)+(cl(4,:)-cl(3,:))/(range(4)-range(3))*(Re-range(3));
      CD=cd(3,:)+(cd(4,:)-cd(3,:))/(range(4)-range(3))*(Re-range(3));
else
      CL=cl(4,:);
      CD=cd(4,:);
end

end