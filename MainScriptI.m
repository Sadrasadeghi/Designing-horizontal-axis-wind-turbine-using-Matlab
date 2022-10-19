clear all
clc
format long

%INbut Data
i=6;                                                    %Group number
D=26.5+(i-1)*3;                                         %Rotor Diameter
R=D/2;                                                  %Rotor radius
V0=8;                                                   %Wind velocity
Nb=3;                                                   %Number of the blades
lambda=7;                                               %Tip speed ratio
Kv=1.5e-5;                                              %Kinematic Viscocity Value
MMair=28.84;                                            %Molecular Mass of air
Temp_amb=288.15;                                        %Ambient temperature (15°C)
p_amb=1*10^5;                                           %Ambient pressure
Rair=8314/MMair;                                             
ro_air=p_amb/(Rair*Temp_amb);                           %Air density
omega=lambda*V0/R;                                      %Omega
% r_root=R_root/R ecc...
r_root=0.40;
r_hub=0.1;
r_primary=0.75;
r_tip=0.95;
R_root=R*r_root;
R_hub=R*r_hub;
R_primary=R*r_primary;
R_tip=R*r_tip;

%Calculations to find xi and Ri
n=150;
Ri=zeros(1,n+1);
Xi=zeros(1,n+1);

for j=0:n
    Ri(1,j+1)=R_hub+j*(R-R_hub)/n;
    Xi(1,j+1)=omega/V0*Ri(1,j+1);
end

%Computations to find a, a_prime, phi_inf, w_inf and the correction factor
%Fc
ai=zeros(1,n+1);
ai_prime=zeros(1,n+1);
phi_inf=zeros(1,n+1);
w_inf=zeros(1,n+1);
f=zeros(1,n+1);
Fc=zeros(1,n+1);

for k=1:n+1
    
    x=[0.3;0.5];
    tol=10^(-5);
    nmax=100;


    F = @(x) [ (4*x(1)-1)*x(2)-(1-3*x(1));x(1)*(1-x(1))-x(2)*(1+x(2))*Xi(k)^2 ];

    JF = @(x) [4*x(2)+3 , 4*x(1)-1; 1-2*x(1), -(Xi(k)^2)*(1+2*x(2))];

    sol=newtonsys(F,JF,x,tol,nmax);

    ai(k)=sol(1);
    ai_prime(k)=sol(2);
    phi_inf(k)=atand((1-ai(k))/(Xi(k)*(1+ai_prime(k))));
    w_inf(k)=(1-ai(k))*(V0/sind(phi_inf(k)));
    f(k)=(Nb/2)*(R-Ri(k))/(Ri(k)*sind(phi_inf(k)));
    Fc(k)=(2/pi)*acos(exp(-f(k)));
        
end

% % Input Data( Here from some predefined files avaialbe we load our input
%datas)

load S818_10.txt;
load S818_20.txt;
load S818_30.txt;
load S818_40.txt;

load S830_10.txt;
load S830_20.txt;
load S830_30.txt;
load S830_40.txt;

load S832_10.txt;
load S832_20.txt;
load S832_30.txt;
load S832_40.txt;


% Interpoliation part
S818_10=interpolation(S818_10);
S818_20=interpolation(S818_20);
S818_30=interpolation(S818_30);
S818_40=interpolation(S818_40);

S830_10=interpolation(S830_10);
S830_20=interpolation(S830_20);
S830_30=interpolation(S830_30);
S830_40=interpolation(S830_40);

S832_10=interpolation(S832_10);
S832_20=interpolation(S832_20);
S832_30=interpolation(S832_30);
S832_40=interpolation(S832_40);

%Optimal Chord Calculation

REi=ones(1,length(Ri))*10000;
cd=zeros(1,length(Ri));
cl=zeros(1,length(Ri));
eff=zeros(1,length(Ri));
betta_inf=zeros(1,length(Ri));
beta_c=zeros(1,length(Ri));
sigma=zeros(1,length(Ri));
c=zeros(1,length(Ri));

for i=1:length(Ri)
    ai_new=ai(i);
    aip_new=ai_prime(i);
    Re=2500000;
    j=0;
    
     while abs(REi(i)-Re)>10000
         
        Re=REi(i);
        j=j+1;
        
        phi_inf(i)=atand((1-ai_new)/(Xi(i)*(1+aip_new)));
        w_inf(i)=(1-ai_new)*V0/sind(phi_inf(i));
        f(i)=(Nb/2)*(R-Ri(i))/(Ri(i)*sind(phi_inf(i)));
        Fc(i)=(2/pi)*acos(exp(-f(i)));

        if Ri(i)<=R_root
            
            [clS818,cdS818,alphaS818]=interp_Re(S818_10,S818_20,S818_30,S818_40,Re);
            [eff(i),p]=max(clS818./cdS818);    % here we find the max eff and its position in the matrix
            cl(i)=clS818(p);
            cd(i)=cdS818(p);
            eff(i)=cl(i)/cd(i);
            betta_inf(i)=alphaS818(p);

            
        elseif  Ri(i) > R_root && Ri(i) <= R_primary 
            
            [clS818,cdS818,alphaS818]=interp_Re(S818_10,S818_20,S818_30,S818_40,Re);
            [clS830,cdS830,alphaS830]=interp_Re(S830_10,S830_20,S830_30,S830_40,Re);
            [cl(i),cd(i),betta_inf(i)]=interp_Radius(clS818,cdS818,clS830,cdS830,alphaS818,alphaS830,Ri(i),R_root,R_primary);
            eff(i)=cl(i)/cd(i);

      
             elseif Ri(i) > R_primary && Ri(i) <= R_tip
                 
                [clS830,cdS830,alphaS830]=interp_Re(S830_10,S830_20,S830_30,S830_40,Re);
                [clS832,cdS832,alphaS832]=interp_Re(S832_10,S832_20,S832_30,S832_40,Re);
                [cl(i),cd(i),betta_inf(i)]=interp_Radius(clS830,cdS830,clS832,cdS832,alphaS830,alphaS832,Ri(i),R_primary,R_tip);
                eff(i)=cl(i)/cd(i);

            
        else
            [clS832,cdS832,alphaS832]=interp_Re(S832_10,S832_20,S832_30,S832_40,Re);
            [eff(i),p]=max(clS832./cdS832);
            cl(i)=clS832(p);
            cd(i)=cdS832(p);
            eff(i)=cl(i)/cd(i);
            betta_inf(i)=alphaS832(p);

        end
        
    beta_c(i)=phi_inf(i)-betta_inf(i);    
    sigma(i)=4*V0*Fc(i)*ai_new/(w_inf(i)*(cl(i)/tand(phi_inf(i))+cd(i)));
    aip_new=sigma(i)*w_inf(i)/(4*omega*Ri(i)*Fc(i))*(cl(i)-cd(i)/tand(phi_inf(i)));
    c(i)=sigma(i)*2*pi*Ri(i)/Nb;
    ai_new=(1+aip_new)/(4*aip_new+3);
    REi(i)=w_inf(i)*c(i)/Kv;
    
     end
end





chord_corr=zeros(length(c),1);  
sigma_corr=zeros(1,length(c));
ai_final=zeros(1,length(c));
aip_final=zeros(1,length(c));
beta_c_final=zeros(1,length(c));
ai_bem=zeros(1,length(c));
aip_bem=zeros(1,length(c));

coefficients=polyfit(Ri(1,20:150),c(1,20:150),2);
retta= @(x) coefficients(1)*x^2+coefficients(2)*x+coefficients(3);

for k=1:length(c)
    
    chord_corr(k)=retta(Ri(k));
   
    sigma_corr(k)=c(k)*Nb/(2*pi*Ri(k));
    ai_final(k)=sigma(k)*w_inf(k)/(4*V0*Fc(k))*(cl(k)/tand(phi_inf(k))+cd(k));
    aip_final(k)=sigma(k)*w_inf(k)/(4*omega*Ri(k)*Fc(k))*(cl(k)-cd(k)/tand(phi_inf(k)));
end

figure()
plot(Ri/R,c/R);
xlabel('r/R')
ylabel('c/R')
hold on;
plot(Ri/R,chord_corr/R);
hold off;

figure()
plot(Ri,c);
axis equal;

phi_final=zeros(1,length(c));
Winf_final=zeros(1,length(c));
Fc_final=zeros(1,length(c));
f_final=zeros(1,length(c));
Re_final=zeros(1,length(c));
toller=ones(1,length(c))*0.0001;
counter=zeros(1,length(k));

for k=1:length(c)
    
    pippo=ai_final(k);
    paperino=aip_final(k);
    
    counter(k)=0;
    
    while (abs(ai_final(k)-ai_bem(k))>toller(k))||(abs(aip_final(k)-aip_bem(k))>toller(k))
        
        counter(k)=counter(k)+1;
        phi_final(k)=atand((1-pippo)/((1+paperino)*Xi(k)));
        Winf_final(k)=(1-pippo)*V0/sind(phi_final(k));
        f_final(k)=(Nb/2)*(R-Ri(k))/(Ri(k)*sind(phi_final(k)));
        Fc_final(k)=(2/pi)*acos(exp(-f_final(k)));
        Re_final(k)=Winf_final(k)*chord_corr(k)/Kv;
    
        if Ri(k)<=R_root

                [clS818,cdS818,alphaS818]=interp_Re(S818_10,S818_20,S818_30,S818_40,Re_final(k));
                [eff(k),p]=max(clS818./cdS818);    % here we find the max eff and its position in the matrix
                cl(k)=clS818(p);
                cd(k)=cdS818(p);
                eff(k)=cl(k)/cd(k);
                betta_inf(k)=alphaS818(p);


            elseif  Ri(k) > R_root && Ri(k) <= R_primary 

                [clS818,cdS818,alphaS818]=interp_Re(S818_10,S818_20,S818_30,S818_40,Re_final(k));
                [clS830,cdS830,alphaS830]=interp_Re(S830_10,S830_20,S830_30,S830_40,Re_final(k));
                [cl(k),cd(k),betta_inf(k)]=interp_Radius(clS818,cdS818,clS830,cdS830,alphaS818,alphaS830,Ri(k),R_root,R_primary);
                eff(k)=cl(k)/cd(k);


                 elseif Ri(k) > R_primary && Ri(k) <= R_tip

                    [clS830,cdS830,alphaS830]=interp_Re(S830_10,S830_20,S830_30,S830_40,Re_final(k));
                    [clS832,cdS832,alphaS832]=interp_Re(S832_10,S832_20,S832_30,S832_40,Re_final(k));
                    [cl(k),cd(k),betta_inf(k)]=interp_Radius(clS830,cdS830,clS832,cdS832,alphaS830,alphaS832,Ri(k),R_primary,R_tip);
                    eff(k)=cl(k)/cd(k);


            else
                [clS832,cdS832,alphaS832]=interp_Re(S832_10,S832_20,S832_30,S832_40,Re_final(k));
                [eff(k),p]=max(clS832./cdS832);
                cl(k)=clS832(p);
                cd(k)=cdS832(p);
                eff(k)=cl(k)/cd(k);
                betta_inf(k)=alphaS832(p);

        end
       beta_c_final(k)=phi_final(k)-betta_inf(k);
       ai_bem(k)=sigma_corr(k)*Winf_final(k)/(4*V0*Fc_final(k))*(cl(k)/tand(phi_final(k))+cd(k));
       aip_bem(k)=sigma_corr(k)*Winf_final(k)/(4*omega*Ri(k)*Fc_final(k))*(cl(k)-cd(k)/tand(phi_final(k)));
       
       ai_final(k)=pippo;
       aip_final(k)=paperino;
       
       eps=0.1;
       paperino=eps*aip_bem(k)+(1-eps)*aip_final(k);
       pippo=eps*ai_bem(k)+(1-eps)*ai_final(k);
       
    end
end

Integral1=0;
for k=1:(length(c)-1)
    f1= @(r) c(k)*((Winf_final(k))^2)*(cl(k)*sind(phi_final(k))-cd(k)*cosd(phi_final(k)))*r;
    I1 = simpcomp(Ri(k), Ri(k+1), 100, f1);
    Integral1=Integral1+I1;
end
Cp_torque=((Nb*omega)/(pi*(R^2)*(V0^3)))*Integral1;
W_torque=ro_air*Nb*omega/2*Integral1;

Integral2=0;
for k=1:(length(c)-1)
    f2= @(y) Fc_final(k)*aip_bem(k)*(1-ai_bem(k))*y.^3;
    I2 = simpcomp(Xi(k), Xi(k+1), 100, f2);
    Integral2=Integral2+I2;
end
Cp_momentum=(8/(lambda^2))*Integral2;
W_momentum=4*ro_air*pi*(R^2)*(V0^3)/(lambda^2)*Integral2;

figure()
plot(Ri/R,cl);
xlabel('r/R')
ylabel('cl')

figure()
plot(Ri/R,cd);
xlabel('r/R')
ylabel('cd')

figure()
plot(Ri/R,phi_final);
xlabel('r/R')
ylabel('\phi_{final}')

figure()
plot(Ri/R,sigma_corr.*cl);
xlabel('r/R')
ylabel('\sigma_{corr}*cl')

figure()
plot(Ri/R,beta_c_final);
xlabel('r/R')
ylabel('\beta_{cfinal}')

figure()
plot(Xi,ai_bem);
xlabel('x')
ylabel('a_{end},a\prime_{end}'),
hold on;
plot(Xi,aip_bem);

figure()
L=zeros(length(Ri),3);
for i=1:length(Ri)
   L(i,1)=chord_corr(i)*cosd(-beta_c_final(i));
   L(i,2)=chord_corr(i)*sind(-beta_c_final(i));
   L(i,3)=Ri(i);
   plot3([0 L(i,1)],[0 L(i,2)],[1 1]*L(i,3),'k')
   title('Twisted blade')
   axis equal
   hold on
end