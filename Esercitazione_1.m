%Helicopter Flight Mechanics
%Esercitazione 1A - Lynx
%----------- Created by Matheus Padilha -----------%
clear all
close all
clc
%----------- Input -----------%
g=9.81; %[m/s^2] gravity @sl
rho_sl=1.225; %[kg/m^3] density @sl
nb=4; %blade number
c=0.391; %[m] mean aerodynamic chord
R=6.4; %[m] blade radius
Omega=35.63; %[rad/s] angular velocity - main rotor
W_TO=4313.7*g; %[N] operative mass
ne=2; %number of engines
eta=1;%we assume power efficiency to be 100%
Pa_sl=850; %[bhp] engine power @sl in British Horse Power -> must be converted
bhp_w=745.7; %conversion factor between bhp and W
k=1.16; %induced power factor 
Cd0 = 0.008; %profile drag coefficient
R_tr=1.106; %[m] tail rotor radius
x_tr=7.660; %[m] distance from main rotor to tail rotor
Cd0_tr=Cd0;
sigma_tr=0.208;%solidity ratio tail rotor - table
%c_tr=R_tr/7; %[m] mean aerodynamic chord - tail rotor

%------------------------------
%% PART 1
%We suppose the initial weight W_TO

T_SL=288.15; %[K] Temperature at sea level
a=-6.5/1000; %[K/m] temperature gradient
h=linspace(0,10000,11); %[m] height distribution - 6 points
T=T_SL+a.*h; %[K] temperature profile 

rho=rho_sl*(T./T_SL).^4.25; %[kg/m^3] air density with altitude
Pa=ne*Pa_sl*bhp_w*eta.*rho./rho_sl.*(1/1000); %[kW] available power for 2 engines

%-------- Display -------%
disp('Helicopter Flight Mechanics')
disp('Esercitazione 1 - Lynx')
disp(' ')
disp('----Part 1 - Results for variable altitude from 0 to 10000 m----')
disp(' ')
disp('Weight (operative weight) [N]:')
disp(W_TO)
disp('Density [kg/m^3]:')
disp(rho)
disp('Available power Pd [kW]:')
disp(Pa) %this way we may represent in output the power in kW

figure()
plot(h,Pa,'LineWidth',2)
xlabel('Altitude [m]')
ylabel('Available power [kW]')
grid on

%%
%Main rotor power
%Induced power
A=pi*R^2;
sigma=(nb*c*R)/A;%solidity ratio
v_star=sqrt(W_TO/A).*sqrt(1./(2*rho)); %[m/s] induced velocity
Pi_star=(k*W_TO.*v_star)./1000; %[kW] effective induced power in hover


%%so far using TDA we have ignored the power due to other losses such as
%%the profile power Po
Po_star=(rho.*(Omega*R)^3*sigma*A*Cd0)./(8*1000); %[kW] profile power in hover
P_mr=Pi_star+Po_star; %[kW] Main rotor's power in hover
FM=Pi_star/(Pi_star+Po_star); %figure of merit

%%
%%Tail rotor power

T_tr=(P_mr*R)./(Omega*R*x_tr)*1000; %[N] thrust needed on tail rotor to equilibrate the torque provided by the main rotor
A_tr=pi*R_tr^2; %[m^2] tail rotor's area


Omega_tr=6*Omega; %[rad/s] tail rotor's angular velocity
Po_tr_star=rho.*A_tr*(Omega_tr*R_tr)^3*(sigma_tr*Cd0_tr)/(8*1000); %[kW] tail rotor's profile power in hover
Pi_tr_star=k*T_tr.*sqrt(T_tr./A_tr).*sqrt(1./(2*rho))*1/1000; %[kW] tail rotor's induced power in hover
P_tr=Po_tr_star+Pi_tr_star; %[kW] tail rotor's power

%%Total power
P_tot=P_mr+P_tr;

dP=Pa-P_tot;
figure()
plot(h,Pa,h,P_tot,h,P_mr,h,P_tr,h,dP, 'LineWidth',2)
xlabel('Altitude [m]')
ylabel('Power [kW]')
grid on
legend('Available power','Total power','Main rotor power','Tail rotor power','Power excess','Location','best')
title ('Power distribution in hover mode')

disp(' ')
disp('Altitude [m]:')
disp(h)
disp('Main rotor profile power - hover [kW]:')
disp(Po_star)
disp('Main rotor induced power - hover [kW]:')
disp(Pi_star)
disp('Tail rotor profile power - hover [kW]:')
disp(Po_tr_star)
disp('Tail rotor induced power - hover [kW]:')
disp(Pi_tr_star)
%% Ceiling (quota di tangenza)
%to find the ceiling condition we shall refine the density profile
h=0:100:10000;

T=T_SL+a.*h; %[K] temperature profile 
rho=rho_sl*(T./T_SL).^4.25; %[kg/m^3] air density with altitude
Pa=ne*Pa_sl*bhp_w*eta.*rho./rho_sl; %[W] available power for 2 engines

%main rotor power
v_star=sqrt(W_TO/A).*sqrt(1./(2*rho)); %[m/s] induced velocity
Pi_star=(k*W_TO.*v_star); %[W] effective induced power
Po_star=(rho.*(Omega*R)^3*sigma*A*Cd0)./(8);%[W] profile power
P_mr=Pi_star+Po_star; %[W] Main rotor's power


%tail rotor power
T_tr=(P_mr*R)./(Omega*R*x_tr); %[N] thrust needed on tail rotor to equilibrate the torque provided by the main rotor
Po_tr_star=rho.*A_tr*(Omega_tr*R_tr)^3*(sigma_tr*Cd0_tr)/(8); %[W] tail rotor's profile power
Pi_tr_star=k*T_tr.*sqrt(T_tr./A_tr).*sqrt(1./(2*rho)); %[W] tail rotor's induced power
P_tr=Po_tr_star+Pi_tr_star; %[W] tail rotor's power

P_tot=P_mr+P_tr;
%to find ceiling I assume the point of inversion the one that has an excess
%power of at most 1 kW
dP=Pa-P_tot;

v=find(abs(dP)<1000,1,'first'); %finds the first value where the excess of power in magnitude is lower than 1kW (positive)
h_c=h(v); %[m] ceiling quote
disp('Ceiling [m]:')
disp(h_c)

figure()
hold on 
plot(h,P_tot/1000,h,Pa/1000,h,Po_star/1000,h,Pi_star/1000,h,P_tr/1000, 'LineWidth',2)
xlabel('Altitude [m]')
ylabel('Power [W]')
grid on
plot(h(v),P_tot(v)/1000,'*r')
legend('Total power','Available power','Profile power - main rotor','Induced power - main rotor','Tail rotor power','ceiling','Location','best')
title ('Power distribution and ceiling condition - hover ')
hold off

%% Part 2
%evaluate the required power in forward flight @sl 
rho=rho_sl; %[kg/m^3] density
K=4.7;% forward flight coefficient
V=0:1:120;%[m/s] forward flight velocities
Pa=(ne*Pa_sl*bhp_w)/(1000)*ones(1,length(V)); %[kW] available power = available power @sl

mu=V./(Omega*R); %rapporto di avanzamento 

Cp0_1=(sigma*Cd0)/8*(1+K.*mu.^2); %profile drag coefficient for low forward velocity (mu<0.2)
Cp0_2=(sigma*Cd0)/8*(1+4*mu.^2+5/8*mu.^4);%profile drag coefficient for high forward velocity (mu>0.2)

%Profile power (main rotor)
Po_1=(rho*A*(Omega*R)^3.*Cp0_1)./(1000); %[kW]profile power with no reverse flow effect
Po_2=(rho*A*(Omega*R)^3.*Cp0_2)./(1000);  %[kW]profile power with reverse flow effect
% for i=1:length(mu)
%     if mu(i)<0.2
%         Po_2(i)=(rho*A*(Omega*R)^3.*Cp0_1(i))/1000;%[kW]profile power with no reverse flow effect as mu<0.2
%     else
%         Po_2(i)=(rho*A*(Omega*R)^3.*Cp0_2(i))/1000;%[kW]profile power with reverse flow effect as mu>=0.2
%     end
% end

% figure()
% plot(mu,Po_1,mu,Po_2)
% legend('1','2')

%Induced power (main rotor)
%We suppose T=W_TO
v_star=sqrt(W_TO/A).*sqrt(1./(2*rho)); %[m/s] induced velocity in hover
vi=sqrt(sqrt(V.^4./4+v_star^4)-V.^2./2); %[m/s] induced velocity in forward flight

Pi=(W_TO.*vi)./(1000); %[kW] induced power for forward flight

%Parassite power
kf=1.2;
cf=0.0040;
kw=2.04;
f=kf*cf*kw*A;
P_p=1/2*rho*f*V.^3*1e-3; %[kW] parassite power

%Tail rotor power
T_tr=(k*W_TO)/(Omega*x_tr)*sqrt(W_TO/(2*rho*A)); %[N] Tail rotor thrust
DL_tr=T_tr/A_tr; %[N/m^2] disk load for the tail rotor

u0_tr=sqrt(DL_tr/(2*rho)); %[m/s] induced velocity on tail rotor

P_tr=T_tr*u0_tr*1e-3*ones(1,length(V)); %[kW] Tail rotor's (induced) power on forward flight
%In forward flight we do not compute the profile power for the tail rotor

%Other powers
%By convention we consider the other powers due to: losses due to non
%uniform velocty profile, tip losses and wake downstream the rotor
Pi_star=k*W_TO*sqrt(W_TO/(2*rho*A))*1e-3; %[kW] induced power (main rotor) in hover
P_other=0.17*Pi_star*ones(1,length(V));

P_mr=Po_1+Pi+P_other; %[kW] main rotor power

P_tot=P_mr+P_p+P_tr; %[kW] total power
figure()
plot(V,P_tr,V,Pi,V,Po_1,V,Po_2,V,P_p,V,P_other,V,P_tot,V,Pa,'--','LineWidth',1.5)
legend('Tail rotor power','Induced power','Profile power - no RF','Profile power - with RF','Parassite power','Other power','Total power','Available power @sl','Location','best')
grid on
xlabel ('Forward velocity V [m/s]')
ylabel ('Power [kW]')
title ('Power distribution - forward flight')

% figure()
% plot(V,P_tr,V,P_mr,V,P_tot,V,Pa,'LineWidth',2)
% legend('Tail rotor power','Main rotor power - no RF','Total power','Available power @sl','Location','best')
% grid on
% xlabel ('Forward velocity V [m/s]')
% ylabel ('Power [kW]')
% title ('Power distribution - forward flight')

% Reverse flow error
error=(abs(Cp0_2-Cp0_1)./Cp0_2).*100; %[%] error for not considering the reverse flow;
figure ()
plot(mu,error,'LineWidth',1.5)
xlabel('Non-dimensional forward flight velocity \mu')
ylabel('Cp_0 Relative error [%]')
set(gca,'Xlim',[0,mu(end)])
grid on
%%
%plot difference in power when considering reverse flow or not
P_mr_1=Po_1+Pi+P_other; %[kW] main rotor power - no RF
P_mr_2=Po_2+Pi+P_other; %[kW] main rotor power - w/ RF
P_tot_1=P_mr_1+P_p+P_tr; %[kW] total power - no RF
P_tot_2=P_mr_2+P_p+P_tr; %[kW] total power - w/ RF

figure()
plot(V,P_tot_1,'k',V,P_tot_2,'r','LineWidth',1.5)
title('Total power - Reverse flow effect')
xlabel('Forward velocity V [m/s')
ylabel('Pwer [kW]')
legend ('No reverse flow', 'With reverse flow')