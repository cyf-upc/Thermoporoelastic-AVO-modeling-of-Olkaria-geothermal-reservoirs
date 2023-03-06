clear;clc;
%%
f=1.0*power(10,9);%频率
omega=2*pi*f;
sigma=1; 
theta=1; 
gg=0; gamma0=0/180*pi; %gg=0,LS-theory,实线；gg=1，GL-theory，虚线%衰减角
% gg=0; gamma0=40/180*pi;  %gg=0,LS-theory,实线；gg=1，GL-theory，虚线%衰减角
% gg=1; gamma0=0/180*pi; %gg=0,LS-theory,实线；gg=1，GL-theory，虚线%衰减角
% gg=1; gamma0=40/180*pi; %gg=0,LS-theory,实线；gg=1，GL-theory，虚线%衰减角
ZI=0;%
for ii=1:1:2
    if ii==1
%% 入射地层
Ks(ii)=35*10^9;  rhos(ii)=2750;      Km(ii)=10*10^9;   mu(ii)=12*10^9;  
phi(ii)=0.05;     bar_kappa(ii)=1*9.87*10^(-12);         tort(ii)=2;       
rhof(ii)=1000;    eta(ii)=0.001;       Kf(ii)=2.4*10^9; betaf(ii)=0.8*10^6; 
c(ii)=3.2*10^5;  beta(ii)=1.2*10^5;  T0(ii)=600;        bar_gamma(ii)=1;

tau1(ii)=0.3*10^-9;  tau2(ii)=0.3*10^-9; tau3(ii)=0.15*10^-9;  tau4(ii)=0; %GL模型
tau1(ii)=0;          tau2(ii)=0;         tau3(ii)=0.15*10^-9;  tau4(ii)=0.15*10^-9;      %LS模型
    elseif ii==2
%% 透射地层
Ks(ii)=45*10^9;  rhos(ii)=2850;      Km(ii)=11*10^9;   mu(ii)=13*10^9;  
phi(ii)=0.1;     bar_kappa(ii)=1*9.87*10^(-12);         tort(ii)=2;       
rhof(ii)=1000;   eta(ii)=0.001;      Kf(ii)=2.4*10^9;   betaf(ii)=0.8*10^6; 
c(ii)=8.2*10^5;  beta(ii)=2.4*10^5;  T0(ii)=650;        bar_gamma(ii)=2.6;

tau1(ii)=0.3*10^-9;  tau2(ii)=0.3*10^-9; tau3(ii)=0.15*10^-9;  tau4(ii)=0; %GL模型
tau1(ii)=0;  tau2(ii)=0;  tau3(ii)=0.15*10^-9; tau4(ii)=0.15*10^-9;      %LS模型
    end
    
lambda(ii)=Km(ii)-2*mu(ii)/3;  
b(ii)=eta(ii)/bar_kappa(ii);      E(ii)=lambda(ii)+2*mu(ii); 
m(ii)=tort(ii)*rhof(ii)/phi(ii);  rho(ii)=(1-phi(ii))*rhos(ii)+phi(ii)*rhof(ii);
MM(ii)=Ks(ii)/(1-phi(ii)-Km(ii)/Ks(ii)+phi(ii)*Ks(ii)/Kf(ii));  
bar_alpha(ii)=1-Km(ii)/Ks(ii);    EG(ii)=E(ii)+bar_alpha(ii)^2*MM(ii); 

bar_tau1(ii)=1-1i*tau1(ii)*omega;  bar_tau2(ii)=1-1i*tau2(ii)*omega; 
bar_tau3(ii)=1-1i*tau3(ii)*omega;  bar_tau4(ii)=1-1i*tau4(ii)*omega;

coe6(ii) =1i*omega^2*phi(ii)*bar_gamma(ii)*MM(ii)*E(ii);  
coe41(ii)=b(ii)*bar_gamma(ii)*EG(ii)-1i*omega*bar_gamma(ii)*(m(ii)*EG(ii)+MM(ii)*(rho(ii)-2*bar_alpha(ii)*rhof(ii)))+c(ii)*MM(ii)*E(ii)*bar_tau3(ii); %
coe42(ii)=betaf(ii)*E(ii)*bar_tau2(ii)+(1-bar_alpha(ii))*MM(ii)*(phi(ii)*beta(ii)*bar_tau1(ii)-bar_alpha(ii)*betaf(ii)*bar_tau2(ii)); %
coe4(ii) =omega^2*(omega*(phi(ii)*coe41(ii)+beta(ii)*T0(ii)*bar_tau4(ii)*coe42(ii))); 
coe21(ii)=omega*rho(ii)*bar_gamma(ii)+1i*c(ii)*EG(ii)*bar_tau3(ii)+1i*beta(ii)^2*T0(ii)*bar_tau1(ii)*bar_tau4(ii);%(w(βΦt1(m-ρf)+t2βf(ρ-ρf)))
coe22(ii)=m(ii)*EG(ii)+MM(ii)*(rho(ii)-2*bar_alpha(ii)*rhof(ii));
coe23(ii)=beta(ii)*phi(ii)*(m(ii)-rhof(ii))*bar_tau1(ii)+betaf(ii)*(rho(ii)-rhof(ii))*bar_tau2(ii); %(Em+Ma^2m+M(rho-2aρf))
coe2(ii) =omega^4*(-b(ii)*phi(ii)*coe21(ii)+omega*(1i*omega*bar_gamma(ii)*phi(ii)*(m(ii)*rho(ii)-rhof(ii)^2)-c(ii)*phi(ii)*bar_tau3(ii)*coe22(ii)-beta(ii)*T0(ii)*bar_tau4(ii)*coe23(ii))); %波数二次方系数 c
coe0(ii) =omega^6*(c(ii)*phi(ii)*bar_tau3(ii)*(1i.*b(ii)*rho(ii)-omega*(rhof(ii)^2-m(ii)*rho(ii)))); %d

coe(ii,:)=[coe6(ii) coe4(ii) coe2(ii) coe0(ii)]; 
k_ptb(ii,:)=sort(sqrt(roots(coe(ii,:))));

k_S(ii)=sqrt(omega^2.*(1i.*b(ii).*rho(ii)-m(ii).*rho(ii).*omega+omega.*rhof(ii).^2)./(mu(ii).*(1i.*b(ii)-omega.*m(ii))));%横波波数

end

k(1)=k_ptb(1,1); k(2)=k_ptb(1,2); k(3)=k_ptb(1,3); k(4)=k_S(1);  
k(5)=k_ptb(2,1); k(6)=k_ptb(2,2); k(7)=k_ptb(2,3); k(8)=k_S(2);

k0=sigma*k(1)+(1-sigma)*k(4);

%%
MN=1000;
dg=90/MN;

for ii=1:MN
 
 alfa=(ii-1)*dg/180*pi;  ppp(ii)=alfa/pi*180;   
  
 P0=sqrt(0.5*( real(k0^2)+sqrt((real(k0^2))^2+(imag(k0^2)/cos(gamma0))^2))); 
 A0=sqrt(0.5*(-real(k0^2)+sqrt((real(k0^2))^2+(imag(k0^2)/cos(gamma0))^2)));
 p0=(P0)*sin(alfa)+1i*(A0)*sin(alfa-gamma0);

for j=1:3  
    p(j)=p0;
    q(j)=(sqrt(k(j)^2-p0^2));
    nu(j)=( (omega/k(j))^2*( rho(1)*betaf(1)*bar_tau2(1)-rhof(1)*phi(1)*beta(1)*bar_tau1(1) )+bar_alpha(1)*MM(1)*phi(1)*beta(1)*bar_tau1(1)-EG(1)*betaf(1)*bar_tau2(1) )./...
          ( (omega/k(j))^2*( phi(1)*beta(1)*bar_tau1(1)*(m(1)+1i*b(1)/omega)-betaf(1)*rhof(1)*bar_tau2(1) )+bar_alpha(1)*MM(1)*betaf(1)*bar_tau2(1)-phi(1)*MM(1)*beta(1)*bar_tau1(1) );%中间参数                  
    delta(j)=( T0(1)*beta(1)*bar_tau4(1)*omega^2.*(1+nu(j)) )/( c(1)*bar_tau3(1)*(omega/k(j))^2+1i*bar_gamma(1)*omega );    %中间参数，温度
end

for j=5:7  
    p(j)=p0;
    q(j)=sqrt(k(j)^2-p0^2);
    nu(j)=( (omega/k(j))^2*( rho(2)*betaf(2)*bar_tau2(2)-rhof(2)*phi(2)*beta(2)*bar_tau1(2) )+bar_alpha(2)*MM(2)*phi(2)*beta(2)*bar_tau1(2)-EG(2)*betaf(2)*bar_tau2(2) )./...
          ( (omega/k(j))^2*( phi(2)*beta(2)*bar_tau1(2)*(m(2)+1i*b(2)/omega)-betaf(2)*rhof(2)*bar_tau2(2) )+bar_alpha(2)*MM(2)*betaf(2)*bar_tau2(2)-phi(2)*MM(2)*beta(2)*bar_tau1(2) );%中间参数                  
    delta(j)=( T0(2)*beta(2)*bar_tau4(2)*omega^2.*(1+nu(j)) )/( c(2)*bar_tau3(2)*(omega/k(j))^2+1i*bar_gamma(2)*omega );    %中间参数，温度
end
 nu(4)=-rhof(1)/(m(1)+1i*b(1)/omega);  delta(4)=0; p(4)=p0; 
 q(4)=sqrt(k(4)^2-p0^2);
 nu(8)=-rhof(2)/(m(2)+1i*b(2)/omega);  delta(8)=0; p(8)=p0; 
 q(8)=sqrt(k(8)^2-p0^2);
 q0=sigma*(sqrt(k0^2-p0^2))+(1-sigma)*(sqrt(k(4)^2-p0^2));
 nu0=sigma*nu(1)+(1-sigma)*nu(4);  delta0=sigma*delta(1)+(1-sigma)*delta(4);
               
 %%
              M11=q(1);  M12=q(2);  M13=q(3);  M14=p(4);
              M15=q(5);  M16=q(6);  M17=q(7);  M18=-p(8);
              
              M21=p(1);  M22=p(2);  M23=p(3);  M24=-q(4);
              M25=-p(5); M26=-p(6); M27=-p(7); M28=-q(8);
              
              M31=-( ( lambda(1)+bar_alpha(1)^2*MM(1)+bar_alpha(1)*MM(1)*nu(1) )*k(1)^2+2*mu(1)*q(1)^2+beta(1)*delta(1)*bar_tau1(1) ); 
              M32=-( ( lambda(1)+bar_alpha(1)^2*MM(1)+bar_alpha(1)*MM(1)*nu(2) )*k(2)^2+2*mu(1)*q(2)^2+beta(1)*delta(2)*bar_tau1(1) );  
              M33=-( ( lambda(1)+bar_alpha(1)^2*MM(1)+bar_alpha(1)*MM(1)*nu(3) )*k(3)^2+2*mu(1)*q(3)^2+beta(1)*delta(3)*bar_tau1(1) );  
              M34=-2*mu(1)*p(4)*q(4);  
              M35=( ( lambda(2)+bar_alpha(2)^2*MM(2)+bar_alpha(2)*MM(2)*nu(5) )*k(5)^2+2*mu(2)*q(5)^2+beta(2)*delta(5)*bar_tau1(2));  
              M36=( ( lambda(2)+bar_alpha(2)^2*MM(2)+bar_alpha(2)*MM(2)*nu(6) )*k(6)^2+2*mu(2)*q(6)^2+beta(2)*delta(6)*bar_tau1(2));  
              M37=( ( lambda(2)+bar_alpha(2)^2*MM(2)+bar_alpha(2)*MM(2)*nu(7) )*k(7)^2+2*mu(2)*q(7)^2+beta(2)*delta(7)*bar_tau1(2));   
              M38=-2*mu(2)*p(8)*q(8);
              
              M41=2*mu(1)*p0*q(1); 
              M42=2*mu(1)*p0*q(2);
              M43=2*mu(1)*p0*q(3);
              M44=(p0^2-q(4)^2)*mu(1);
              M45=2*mu(2)*p0*q(5);
              M46=2*mu(2)*p0*q(6);
              M47=2*mu(2)*p0*q(7);
              M48=-(p0^2-q(8)^2)*mu(2);
              
              M51=q(1)*phi(1)*(nu(1)-1); 
              M52=q(2)*phi(1)*(nu(2)-1);
              M53=q(3)*phi(1)*(nu(3)-1);
              M54=p(4)*phi(1)*(nu(4)-1);
              M55=theta*(q(5)*phi(2)*(nu(5)-1));
              M56=theta*(q(6)*phi(2)*(nu(6)-1));
              M57=theta*(q(7)*phi(2)*(nu(7)-1));
              M58=theta*(-p(8)*phi(2)*(nu(8)-1));
              
              M61=theta*( MM(1)*k(1)^2*(bar_alpha(1)+nu(1))+betaf(1)*delta(1)*bar_tau2(1)/phi(1) );  
              M62=theta*( MM(1)*k(2)^2*(bar_alpha(1)+nu(2))+betaf(1)*delta(2)*bar_tau2(1)/phi(1) );
              M63=theta*( MM(1)*k(3)^2*(bar_alpha(1)+nu(3))+betaf(1)*delta(3)*bar_tau2(1)/phi(1) );
              M64=0;
              M65=theta*(-( MM(2)*k(5)^2*(bar_alpha(2)+nu(5))+betaf(2)*delta(5)*bar_tau2(2)/phi(2) )+ZI*omega*q(5)*phi(2)*(nu(5)-1))+(1-theta)*(q(5)*(nu(5)-1));
              M66=theta*(-( MM(2)*k(6)^2*(bar_alpha(2)+nu(6))+betaf(2)*delta(6)*bar_tau2(2)/phi(2) )+ZI*omega*q(6)*phi(2)*(nu(6)-1))+(1-theta)*(q(6)*(nu(6)-1));
              M67=theta*(-( MM(2)*k(7)^2*(bar_alpha(2)+nu(7))+betaf(2)*delta(7)*bar_tau2(2)/phi(2) )+ZI*omega*q(7)*phi(2)*(nu(7)-1))+(1-theta)*(q(7)*(nu(7)-1));
              M68=theta*(-ZI*omega*p(8)*phi(2)*(nu(8)-1))+(1-theta)*(-p(8)*phi(2)*(nu(8)-1));
            
              M71=delta(1); 
              M72=delta(2);
              M73=delta(3);
              M74=0;
              M75=-delta(5);
              M76=-delta(6);
              M77=-delta(7);
              M78=0;
              
              M81=bar_gamma(1)*delta(1)*q(1);  
              M82=bar_gamma(1)*delta(2)*q(2);
              M83=bar_gamma(1)*delta(3)*q(3);
              M84=0;
              M85=bar_gamma(2)*delta(5)*q(5);
              M86=bar_gamma(2)*delta(6)*q(6);
              M87=bar_gamma(2)*delta(7)*q(7);
              M88=0;
              
              N1=q0*sigma-(1-sigma)*p0;
              N2=-p0*sigma-(1-sigma)*q0;
              N3=( ( lambda(1)+bar_alpha(1)^2*MM(1)+bar_alpha(1)*MM(1)*nu0 )*(k0^2)+2*mu(1)*q0^2+beta(1)*delta0*bar_tau1(1) )*sigma...
                  -(1-sigma)*2*mu(1)*p0*q0;  
              N4=2*mu(1)*p0*q0*sigma+(1-sigma)*mu(1)*(q0^2-p0^2);
                  
              N5=q0*phi(1)*(nu0-1)*sigma...
                  -(1-sigma)*p0*phi(1)*(nu0-1); 
              N6=theta*(-( MM(1)*(k0^2)*(bar_alpha(1)+nu0)+betaf(1)*delta0*bar_tau2(1)/phi(1) )*sigma+(1-sigma)*0);
              N7=-delta0*sigma+(1-sigma)*0; 
              N8=bar_gamma(1)*delta0*q0*sigma+(1-sigma)*0; 
              

M_88=[M11,M12,M13,M14,M15,M16,M17,M18;...
      M21,M22,M23,M24,M25,M26,M27,M28;...
      M31,M32,M33,M34,M35,M36,M37,M38;...
      M41,M42,M43,M44,M45,M46,M47,M48;...
      M51,M52,M53,M54,M55,M56,M57,M58;...
      M61,M62,M63,M64,M65,M66,M67,M68;...
      M71,M72,M73,M74,M75,M76,M77,M78;...
      M81,M82,M83,M84,M85,M86,M87,M88];
N_81=[N1;N2;N3;N4;N5;N6;N7;N8];

AA=inv(M_88)*N_81; 
 

 for j=1:8
     RT_ppts(ii,j)       = abs(AA(j)*k(j)/k0); 
     phase_RT_ppts(ii,j) = angle(AA(j)*k(j)/k0);
 end

end
 %%   
figure(1)
subplot(2,4,1) %多图绘制
if gg==0 && gamma0==0/180*pi;
    plot(ppp, RT_ppts(:,1),'k-', 'LineWidth',2.5, 'MarkerSize',2.5); hold on
elseif gg==0 && gamma0==40/180*pi;
    plot(ppp, RT_ppts(:,1),'r-', 'LineWidth',2.5, 'MarkerSize',2.5);hold on
elseif gg==1 && gamma0==0/180*pi;
    plot(ppp, RT_ppts(:,1),'k--', 'LineWidth',2.5, 'MarkerSize',2.5);hold on
elseif gg==1 && gamma0==40/180*pi; 
    plot(ppp, RT_ppts(:,1),'r--', 'LineWidth',2.5, 'MarkerSize',2.5);hold on
end
xlim([0 90]);
set(gca,'ycolor','k');
xlabel('Incidence angle (\circ)');
ylabel('{\itR}_P_{P1}');
hold on;box off;

subplot(2,4,2)
if gg==0 && gamma0==0/180*pi;
    plot(ppp, RT_ppts(:,2),'k-', 'LineWidth',2.5, 'MarkerSize',2.5); hold on
elseif gg==0 && gamma0==40/180*pi;
    plot(ppp, RT_ppts(:,2),'r-', 'LineWidth',2.5, 'MarkerSize',2.5);hold on
elseif gg==1 && gamma0==0/180*pi;
    plot(ppp, RT_ppts(:,2),'k--', 'LineWidth',2.5, 'MarkerSize',2.5);hold on
elseif gg==1 && gamma0==40/180*pi; 
    plot(ppp, RT_ppts(:,2),'r--', 'LineWidth',2.5, 'MarkerSize',2.5);hold on
end
set(gca,'ycolor','k');
ylabel('{\itR}_P_{P2}');xlim([0 90]);
xlabel('Incidence angle (\circ)');
hold on; box off;

subplot(2,4,3)
if gg==0 && gamma0==0/180*pi;
    plot(ppp, RT_ppts(:,3),'k-', 'LineWidth',2.5, 'MarkerSize',2.5); hold on
elseif gg==0 && gamma0==40/180*pi;
    plot(ppp, RT_ppts(:,3),'r-', 'LineWidth',2.5, 'MarkerSize',2.5);hold on
elseif gg==1 && gamma0==0/180*pi;
    plot(ppp, RT_ppts(:,3),'k--', 'LineWidth',2.5, 'MarkerSize',2.5);hold on
elseif gg==1 && gamma0==40/180*pi; 
    plot(ppp, RT_ppts(:,3),'r--', 'LineWidth',2.5, 'MarkerSize',2.5);hold on
end
set(gca,'ycolor','k');
ylabel('{\itR}_P_T');xlim([0 90]);
xlabel('Incidence angle (\circ)');
hold on; box off;

subplot(2,4,4)
if gg==0 && gamma0==0/180*pi;
    plot(ppp, RT_ppts(:,4),'k-', 'LineWidth',2.5, 'MarkerSize',2.5); hold on
elseif gg==0 && gamma0==40/180*pi;
    plot(ppp, RT_ppts(:,4),'r-', 'LineWidth',2.5, 'MarkerSize',2.5);hold on
elseif gg==1 && gamma0==0/180*pi;
    plot(ppp, RT_ppts(:,4),'k--', 'LineWidth',2.5, 'MarkerSize',2.5);hold on
elseif gg==1 && gamma0==40/180*pi; 
    plot(ppp, RT_ppts(:,4),'r--', 'LineWidth',2.5, 'MarkerSize',2.5);hold on
end
set(gca,'ycolor','k');
ylabel('{\itR}_P_S');xlim([0 90]);
xlabel('Incidence angle (\circ)');
hold on; box off;

subplot(2,4,5)
if gg==0 && gamma0==0/180*pi;
    plot(ppp, RT_ppts(:,5),'k-', 'LineWidth',2.5, 'MarkerSize',2.5); hold on
elseif gg==0 && gamma0==40/180*pi;
    plot(ppp, RT_ppts(:,5),'r-', 'LineWidth',2.5, 'MarkerSize',2.5);hold on
elseif gg==1 && gamma0==0/180*pi;
    plot(ppp, RT_ppts(:,5),'k--', 'LineWidth',2.5, 'MarkerSize',2.5);hold on
elseif gg==1 && gamma0==40/180*pi; 
    plot(ppp, RT_ppts(:,5),'r--', 'LineWidth',2.5, 'MarkerSize',2.5);hold on
end
set(gca,'ycolor','k');
ylabel('{\itT}_P_{P1}');xlim([0 90]);
xlabel('Incidence angle (\circ)');
hold on; box off;

subplot(2,4,6)
if gg==0 && gamma0==0/180*pi;
    plot(ppp, RT_ppts(:,6),'k-', 'LineWidth',2.5, 'MarkerSize',2.5); hold on
elseif gg==0 && gamma0==40/180*pi;
    plot(ppp, RT_ppts(:,6),'r-', 'LineWidth',2.5, 'MarkerSize',2.5);hold on
elseif gg==1 && gamma0==0/180*pi;
    plot(ppp, RT_ppts(:,6),'k--', 'LineWidth',2.5, 'MarkerSize',2.5);hold on
elseif gg==1 && gamma0==40/180*pi; 
    plot(ppp, RT_ppts(:,6),'r--', 'LineWidth',2.5, 'MarkerSize',2.5);hold on
end
set(gca,'ycolor','k');
ylabel('{\itT}_P_{P2}');xlim([0 90]);
xlabel('Incidence angle (\circ)');
hold on; box off;

subplot(2,4,7)
if gg==0 && gamma0==0/180*pi;
    plot(ppp, RT_ppts(:,7),'k-', 'LineWidth',2.5, 'MarkerSize',2.5); hold on
elseif gg==0 && gamma0==40/180*pi;
    plot(ppp, RT_ppts(:,7),'r-', 'LineWidth',2.5, 'MarkerSize',2.5);hold on
elseif gg==1 && gamma0==0/180*pi;
    plot(ppp, RT_ppts(:,7),'k--', 'LineWidth',2.5, 'MarkerSize',2.5);hold on
elseif gg==1 && gamma0==40/180*pi; 
    plot(ppp, RT_ppts(:,7),'r--', 'LineWidth',2.5, 'MarkerSize',2.5);hold on
end
set(gca,'ycolor','k');
ylabel('{\itT}_P_T');
xlim([0,90]);
xlabel('Incidence angle (\circ)');
hold on; box off;

subplot(2,4,8)
if gg==0 && gamma0==0/180*pi;
    plot(ppp, RT_ppts(:,8),'k-', 'LineWidth',2.5, 'MarkerSize',2.5); hold on
elseif gg==0 && gamma0==40/180*pi;
    plot(ppp, RT_ppts(:,8),'r-', 'LineWidth',2.5, 'MarkerSize',2.5);hold on
elseif gg==1 && gamma0==0/180*pi;
    plot(ppp, RT_ppts(:,8),'k--', 'LineWidth',2.5, 'MarkerSize',2.5);hold on
elseif gg==1 && gamma0==40/180*pi; 
    plot(ppp, RT_ppts(:,8),'r--', 'LineWidth',2.5, 'MarkerSize',2.5);hold on
set(gca,'Position',[0.79 0.08 0.19 0.37]) ;%左、下、右、上
end
xlim([0 90]);
set(gca,'ycolor','k');
xlabel('Incidence angle (\circ)');
ylabel('{\itT}_P_S');
hold on; box off;

