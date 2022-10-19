function [CL,CD,alpha]=interpReynolds(t10,t20,t30,t40,Re)
% Reynolds Interpolation
% interpola i valori di CL e CD per numeri di reynolds diversi rispettto a
% quelli nelle tabelle

ret=10e5:10e5:40e5; % punti fissi di reynolds delle tabelle

alpha=t10(:,1); %angoli di incidenza di lavoro

%creo un'unica matrice di cl e cd a diversi RE
cl=[t10(:,2),t20(:,2),t30(:,2),t40(:,2)]'; % ogni colonna è un angolo
cd=[t10(:,3),t20(:,3),t30(:,3),t40(:,3)]'; %ogni riga è un Re diverso

if Re<ret(1)
      CL=cl(1,:);
      CD=cd(1,:);
elseif Re>=ret(1) && Re<ret(2)
      CL=cl(1,:)+(cl(2,:)-cl(1,:))/(ret(2)-ret(1))*(Re-ret(1));
      CD=cd(1,:)+(cd(2,:)-cd(1,:))/(ret(2)-ret(1))*(Re-ret(1));
elseif Re>=ret(2) && Re<ret(3)
      CL=cl(2,:)+(cl(3,:)-cl(2,:))/(ret(3)-ret(2))*(Re-ret(2));
      CD=cd(2,:)+(cd(3,:)-cd(2,:))/(ret(3)-ret(2))*(Re-ret(2));
elseif Re>=ret(3) && Re<ret(4)
      CL=cl(3,:)+(cl(4,:)-cl(3,:))/(ret(4)-ret(3))*(Re-ret(3));
      CD=cd(3,:)+(cd(4,:)-cd(3,:))/(ret(4)-ret(3))*(Re-ret(3));
else
      CL=cl(4,:);
      CD=cd(4,:);
end

end