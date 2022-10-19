format long;

S818_10=xlsread('S818_10.xlsx');
S818_20=xlsread('S818_20.xlsx');
S818_30=xlsread('S818_30.xlsx');
S818_40=xlsread('S818_40.xlsx');

S830_10=xlsread('S830_10.xlsx');
S830_20=xlsread('S830_20.xlsx');
S830_30=xlsread('S830_30.xlsx');
S830_40=xlsread('S830_40.xlsx');

S832_10=xlsread('S832_10.xlsx');
S832_20=xlsread('S832_20.xlsx');
S832_30=xlsread('S832_30.xlsx');
S832_40=xlsread('S832_40.xlsx');

% Interpoliation part
S818_10=interpolation2(S818_10);
S818_20=interpolation2(S818_20);
S818_30=interpolation2(S818_30);
S818_40=interpolation2(S818_40);

S830_10=interpolation2(S830_10);
S830_20=interpolation2(S830_20);
S830_30=interpolation2(S830_30);
S830_40=interpolation2(S830_40);

S832_10=interpolation2(S832_10);
S832_20=interpolation2(S832_20);
S832_30=interpolation2(S832_30);
S832_40=interpolation2(S832_40);

Vcut_in=3;
Vcut_out=20;
V_rated=1.7*V0;
u_max=70;
omega_max=u_max/R;
lambda_max=omega_max*R/V0;
Vc=omega_max*R/lambda;
chord=chord_corr;
Vi=Vcut_in:0.1:Vcut_out;
ai=ai_bem.*ones(length(Vi),length(Ri));
aip=aip_bem.*ones(length(Vi),length(Ri));
beta_c=beta_c_final;

Vi=Vcut_in:0.1:Vcut_out;
icut_in=find(Vi==Vcut_in);
icut_out=find(Vi==Vcut_out);
i_c=find(Vi==Vc);
i_rated=107;

omega_vector=zeros(1,i_c);
toll=0.001;
Cp_torque=zeros(length(Vi),1);
Cp_momentum=zeros(length(Vi),1);
W_torque=zeros(length(Vi),1);
W_momentum=zeros(length(Vi),1);

for i=1:i_c
    omega_vector(i)=Vi(i)*lambda/R;
    
    for j=1:length(Ri)
        a_new=ai(i,j);
        ap_new=aip(i,j);
        a_for=0;
        ap_for=0;
        
        while ((abs(a_for-ai(i,j))>toll)) && ((abs(ap_for-aip(i,j))>toll))

                phi_inf(i,j)=atand((1-a_new)/(Xi(j)*(1+ap_new)));
                W_inf(i,j)=(1-a_new)*Vi(i)/sind(phi_inf(i,j));
                f(j)=(Nb/2)*(R-Ri(j))/(Ri(j)*sind(phi_inf(i,j)));
                Fc(j)=(2/pi)*acos(exp(-f(j)));
                beta_inf(i,j)=phi_inf(i,j)-beta_c(j);
                REi(i,j)=W_inf(i,j)*chord(j)/Kv;

                if Ri(j)<=R_root

                    [clS818,cdS818,alphaS818]=interp_Re(S818_10,S818_20,S818_30,S818_40,REi(i,j));
                    [p]=search_angle(alphaS818,beta_inf(i,j));    
                    cl(j)=clS818(p);
                    cd(j)=cdS818(p);
                    eff(j)=cl(j)/cd(j);



                elseif  Ri(j) > R_root && Ri(j) <= R_primary 

                    [clS818,cdS818,alphaS818]=interp_Re(S818_10,S818_20,S818_30,S818_40,REi(i,j));
                    [clS830,cdS830,alphaS830]=interp_Re(S830_10,S830_20,S830_30,S830_40,REi(i,j));
                    [p1]=search_angle(alphaS818,beta_inf(i,j)); 
                    [p2]=search_angle(alphaS830,beta_inf(i,j)); 
                    [cl(j),cd(j)]=interp_Radius2(clS818(p1),cdS818(p1),clS830(p2),cdS830(p2),Ri(j),R_root,R_primary);
                    eff(j)=cl(j)/cd(j);


                 elseif Ri(j) > R_primary && Ri(j) <= R_tip

                    [clS830,cdS830,alphaS830]=interp_Re(S830_10,S830_20,S830_30,S830_40,REi(i,j));
                    [clS832,cdS832,alphaS832]=interp_Re(S832_10,S832_20,S832_30,S832_40,REi(i,j));
                    [p1]=search_angle(alphaS830,beta_inf(i,j)); 
                    [p2]=search_angle(alphaS832,beta_inf(i,j)); 
                    [cl(j),cd(j)]=interp_Radius2(clS830(p1),cdS830(p1),clS832(p2),cdS832(p2),Ri(j),R_primary,R_tip);
                    eff(j)=cl(j)/cd(j);


                    else
                    [clS832,cdS832,alphaS832]=interp_Re(S832_10,S832_20,S832_30,S832_40,REi(i,j));
                    [p]=search_angle(alphaS832,beta_inf(i,j));
                    cd(j)=cdS832(p);
                    cl(j)=clS832(p);
                    eff(j)=cl(j)/cd(j);

                end
            a_for=a_new;
            ap_for=ap_new;
            ai(i,j)=sigma_corr(j)*W_inf(i,j)/(4*Vi(i)*Fc(j))*(cl(j)/tand(phi_inf(i,j))+cd(j));
           aip(i,j)=sigma_corr(j)*W_inf(i,j)/(4*omega_vector(i)*Ri(j)*Fc(j))*(cl(j)-cd(j)/tand(phi_inf(i,j)));
           eps=0.2;
           a_new=ai(i,j)*eps+(1-eps)*a_for;
           ap_new=aip(i,j)*eps+(1-eps)*ap_for;
            end
    end
    
    
        Integral1=0;
    for k=1:(length(c)-1)
        f1= @(r) chord(k)*((W_inf(i,k))^2)*(cl(k)*sind(phi_inf(i,k))-cd(k)*cosd(phi_inf(i,k)))*r;
        I1 = simpcomp(Ri(k), Ri(k+1), 100, f1);
        Integral1=Integral1+I1;
    end
    Cp_torque(i)=((Nb*omega_vector(i))/(pi*(R^2)*(Vi(i)^3)))*Integral1;
    W_torque(i)=ro_air*Nb*omega_vector(i)/2*Integral1;

    Integral2=0;
    for k=1:(length(c)-1)
        f2= @(y) Fc(k)*aip(i,k)*(1-ai(i,k))*y.^3;
        I2 = simpcomp(Xi(k), Xi(k+1), 100, f2);
        Integral2=Integral2+I2;
    end
    Cp_momentum(i)=(8/(lambda^2))*Integral2;
    W_momentum(i)=4*ro_air*pi*(R^2)*(Vi(i)^3)/(lambda^2)*Integral2;
    
end

delta_beta=zeros(1,length(Vi));
toll2=0.001;

for i=(i_c+1):1:i_rated
    
    stalled=0;
    
    for k=1:length(Ri)
    Xi(k)=omega_max*Ri(k)/Vi(i);
    end
    
    while (stalled/length(Ri)<0.2)
        
        for j=1:length(Ri)
            a_new=ai(i,j);
            ap_new=aip(i,j);
            a_for=0;
            ap_for=0;

            while ((abs(a_for-ai(i,j))>toll2)) && ((abs(ap_for-aip(i,j))>toll2))

                phi_inf(i,j)=atand((1-a_new)/(Xi(j)*(1+ap_new)));
                W_inf(i,j)=(1-a_new)*Vi(i)/sind(phi_inf(i,j));
                f(j)=(Nb/2)*(R-Ri(j))/(Ri(j)*sind(phi_inf(i,j)));
                Fc(j)=(2/pi)*acos(exp(-f(j)));
                beta_inf(i,j)=phi_inf(i,j)-beta_c(j)-delta_beta(i);
                REi(i,j)=W_inf(i,j)*chord(j)/Kv;

                if Ri(j)<=R_root

                    [clS818,cdS818,alphaS818]=interp_Re(S818_10,S818_20,S818_30,S818_40,REi(i,j));
                    [p]=search_angle(alphaS818,beta_inf(i,j));    
                    cl(j)=clS818(p);
                    cd(j)=cdS818(p);
                    eff(j)=cl(j)/cd(j);



                elseif  Ri(j) > R_root && Ri(j) <= R_primary 

                    [clS818,cdS818,alphaS818]=interp_Re(S818_10,S818_20,S818_30,S818_40,REi(i,j));
                    [clS830,cdS830,alphaS830]=interp_Re(S830_10,S830_20,S830_30,S830_40,REi(i,j));
                    [p1]=search_angle(alphaS818,beta_inf(i,j)); 
                    [p2]=search_angle(alphaS830,beta_inf(i,j)); 
                    [cl(j),cd(j)]=interp_Radius2(clS818(p1),cdS818(p1),clS830(p2),cdS830(p2),Ri(j),R_root,R_primary);
                    eff(j)=cl(j)/cd(j);


                 elseif Ri(j) > R_primary && Ri(j) <= R_tip

                    [clS830,cdS830,alphaS830]=interp_Re(S830_10,S830_20,S830_30,S830_40,REi(i,j));
                    [clS832,cdS832,alphaS832]=interp_Re(S832_10,S832_20,S832_30,S832_40,REi(i,j));
                    [p1]=search_angle(alphaS830,beta_inf(i,j)); 
                    [p2]=search_angle(alphaS832,beta_inf(i,j)); 
                    [cl(j),cd(j)]=interp_Radius2(clS830(p1),cdS830(p1),clS832(p2),cdS832(p2),Ri(j),R_primary,R_tip);
                    eff(j)=cl(j)/cd(j);


                    else
                    [clS832,cdS832,alphaS832]=interp_Re(S832_10,S832_20,S832_30,S832_40,REi(i,j));
                    [p]=search_angle(alphaS832,beta_inf(i,j));
                    cd(j)=cdS832(p);
                    cl(j)=clS832(p);
                    eff(j)=cl(j)/cd(j);

                end
            a_for=a_new;
            ap_for=ap_new;
            ai(i,j)=sigma_corr(j)*W_inf(i,j)/(4*Vi(i)*Fc(j))*(cl(j)/tand(phi_inf(i,j))+cd(j));
           aip(i,j)=sigma_corr(j)*W_inf(i,j)/(4*omega_max*Ri(j)*Fc(j))*(cl(j)-cd(j)/tand(phi_inf(i,j)));
           eps=0.2;
           a_new=ai(i,j)*eps+(1-eps)*a_for;
           ap_new=aip(i,j)*eps+(1-eps)*ap_for;
            end
        end

        eff_max=max(eff);
        stalled=0;

        for k=1:length(Ri)
            if eff(k)<(0.95*eff_max)
                stalled=stalled+1;
            end
        end

        delta_beta(i)=delta_beta(i)+0.01;

    end
    
    delta_beta(i)=delta_beta(i)-0.01;
        
    Integral1=0;
        for k=1:(length(c)-1)
            f1= @(r) chord(k)*((W_inf(i,k))^2)*(cl(k)*sind(phi_inf(i,k))-cd(k)*cosd(phi_inf(i,k)))*r;
            I1 = simpcomp(Ri(k), Ri(k+1), 100, f1);
            Integral1=Integral1+I1;
        end
    Cp_torque(i)=((Nb*omega_max)/(pi*(R^2)*(Vi(i)^3)))*Integral1;
    W_torque(i)=ro_air*Nb*omega_max/2*Integral1;

    Integral2=0;
        for k=1:(length(c)-1)
            f2= @(y) Fc(k)*aip(i,k)*(1-ai(i,k))*(y.^3);
            I2 = simpcomp(Xi(k), Xi(k+1), 10, f2);
            Integral2=Integral2+I2;
        end
    Cp_momentum(i)=(8/((omega_max*R/Vi(i))^2))*Integral2;
    W_momentum(i)=4*ro_air*pi*(R^2)*(Vi(i)^3)/((omega_max*R/Vi(i))^2)*Integral2;
        
end
 
figure();
plot(Vi,Cp_torque);

figure();
plot(Vi,W_torque);

figure();
plot(Vi,Cp_momentum);

figure();
plot(Vi,W_momentum);

toll3=1000;
toll4=0.001;


for i=(i_rated+1):1:length(Vi)
    
    for k=1:length(Ri)
    Xi(k)=omega_max*Ri(k)/Vi(i);
    end
    
    
    while (abs(W_torque(i)-W_torque(i_rated))>toll3)

        for j=1:length(Ri)
            
            a_new=ai(i,j);
            ap_new=aip(i,j);
            a_for=0;
            ap_for=0;

            while ((abs(a_for-ai(i,j))>toll4)) && ((abs(ap_for-aip(i,j))>toll4))

                phi_inf(i,j)=real(atand((1-a_new)/(Xi(j)*(1+ap_new))));
                W_inf(i,j)=(1-a_new)*Vi(i)/sind(phi_inf(i,j));
                f(j)=(Nb/2)*(R-Ri(j))/(Ri(j)*sind(phi_inf(i,j)));
                Fc(j)=(2/pi)*acos(exp(-f(j)));
                beta_inf(i,j)=phi_inf(i,j)-beta_c(j)-delta_beta(i);
                REi(i,j)=W_inf(i,j)*chord(j)/Kv;

                if Ri(j)<=R_root

                    [clS818,cdS818,alphaS818]=interp_Re(S818_10,S818_20,S818_30,S818_40,REi(i,j));
                    [p]=search_angle(alphaS818,beta_inf(i,j));    
                    cl(j)=clS818(p);
                    cd(j)=cdS818(p);
                    eff(j)=cl(j)/cd(j);



                elseif  Ri(j) > R_root && Ri(j) <= R_primary 

                    [clS818,cdS818,alphaS818]=interp_Re(S818_10,S818_20,S818_30,S818_40,REi(i,j));
                    [clS830,cdS830,alphaS830]=interp_Re(S830_10,S830_20,S830_30,S830_40,REi(i,j));
                    [p1]=search_angle(alphaS818,beta_inf(i,j)); 
                    [p2]=search_angle(alphaS830,beta_inf(i,j)); 
                    [cl(j),cd(j)]=interp_Radius2(clS818(p1),cdS818(p1),clS830(p2),cdS830(p2),Ri(j),R_root,R_primary);
                    eff(j)=cl(j)/cd(j);


                 elseif Ri(j) > R_primary && Ri(j) <= R_tip

                    [clS830,cdS830,alphaS830]=interp_Re(S830_10,S830_20,S830_30,S830_40,REi(i,j));
                    [clS832,cdS832,alphaS832]=interp_Re(S832_10,S832_20,S832_30,S832_40,REi(i,j));
                    [p1]=search_angle(alphaS830,beta_inf(i,j)); 
                    [p2]=search_angle(alphaS832,beta_inf(i,j)); 
                    [cl(j),cd(j)]=interp_Radius2(clS830(p1),cdS830(p1),clS832(p2),cdS832(p2),Ri(j),R_primary,R_tip);
                    eff(j)=cl(j)/cd(j);


                    else
                    [clS832,cdS832,alphaS832]=interp_Re(S832_10,S832_20,S832_30,S832_40,REi(i,j));
                    [p]=search_angle(alphaS832,beta_inf(i,j));
                    cd(j)=cdS832(p);
                    cl(j)=clS832(p);
                    eff(j)=cl(j)/cd(j);

                end
            a_for=a_new;
            ap_for=ap_new;
            ai(i,j)=sigma_corr(j)*W_inf(i,j)/(4*Vi(i)*Fc(j))*(cl(j)/tand(phi_inf(i,j))+cd(j));
           aip(i,j)=sigma_corr(j)*W_inf(i,j)/(4*omega_max*Ri(j)*Fc(j))*(cl(j)-cd(j)/tand(phi_inf(i,j)));
           eps=0.2;
           a_new=ai(i,j)*eps+(1-eps)*a_for;
           ap_new=aip(i,j)*eps+(1-eps)*ap_for;
            end
        end

        Integral1=0;
                for k=1:(length(c)-1)
                    f1= @(r) chord(k)*((W_inf(i,k))^2)*(cl(k)*sind(phi_inf(i,k))-cd(k)*cosd(phi_inf(i,k)))*r;
                    I1 = simpcomp(Ri(k), Ri(k+1), 100, f1);
                    Integral1=Integral1+I1;
                end
            Cp_torque(i)=((Nb*omega_max)/(pi*(R^2)*(Vi(i)^3)))*Integral1;
            W_torque(i)=ro_air*Nb*omega_max/2*Integral1;
            
                Integral2=0;
        for k=1:(length(c)-1)
            f2= @(y) Fc(k)*aip(i,k)*(1-ai(i,k))*y.^3;
            I2 = simpcomp(Xi(k), Xi(k+1), 100, f2);
            Integral2=Integral2+I2;
        end
    Cp_momentum(i)=(8/((omega_max*R/Vi(i))^2))*Integral2;
    W_momentum(i)=4*ro_air*pi*(R^2)*(Vi(i)^3)/((omega_max*R/Vi(i))^2)*Integral2;
    
            delta_beta(i)=delta_beta(i)+0.01;

    end
    
    delta_beta(i)=delta_beta(i)-0.01;    
        
end
 
figure();
plot(Vi,Cp_torque);

figure();
plot(Vi,W_torque);

figure();
plot(Vi,Cp_momentum);

figure();
plot(Vi,W_momentum);