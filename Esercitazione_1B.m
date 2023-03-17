%Helicopter Flight Mechanics
%Esercitazione 1B - Lynx
%----------- Created by Matheus Padilha -----------%
clear
close all
clc
%----------- Input -----------%
g = 9.81; %[m/s^2] gravity @sl
rho_sl = 1.225; %[kg/m^3] density @sl
nb = 4; %blade number
c = 0.391; %[m] mean aerodynamic chord
R = 6.4; %[m] blade radius
Omega = 35.63; %[rad/s] angular velocity - main rotor
W_TO = 4313.7*g; %[N] operative mass
ne = 2; %number of engines
eta = 1;%we assume power efficiency to be 100%
Pa_sl = 850; %[bhp] engine power @sl in British Horse Power -> must be converted
bhp_w = 745.7; %conversion factor between bhp and W
k = 1.16; %induced power factor
Cd0  =  0.008; %profile drag coefficient
R_tr = 1.106; %[m] tail rotor radius
x_tr = 7.660; %[m] distance from main rotor to tail rotor
Cd0_tr = Cd0;
sigma_tr = 0.208;%solidity ratio tail rotor - table
%c_tr = R_tr/7; %[m] mean aerodynamic chord - tail rotor -> decided not to
%use it as we have directly the solidility ratio
Omega_tr = 6*Omega; %[rad/s] tail rotor's angular velocity
A = pi*R^2;
sigma = (nb*c*R)/A;%solidity ratio
K = 4.7;% forward flight coefficient
A_tr = pi*R_tr^2; %[m^2] tail rotor's area

%% Power profile with variable mass and velocity

W = [W_TO:3000:1.4*W_TO]'; %[N] possible weights during TO

V = [0:1:120];%[m/s] forward flight velocities
Pa = (ne*Pa_sl*bhp_w)/(1000)*ones(1,length(V)); %[kW] available power  =  available power @sl
mu = V./(Omega*R); %rapporto di avanzamento

%We shall value the required power without considering the reverse flow
%effect
Cp0 = (sigma*Cd0)/8*(1+K.*mu.^2); %profile drag coefficient without reverse flow (mu<0.2)
Po = (rho_sl*A*(Omega*R)^3.*Cp0)./(1000); %[kW]profile power with no reverse flow effect

v_star = [sqrt(W./A).*sqrt(1./(2*rho_sl))]; %[m/s] induced velocity in hover

for i = 1:length(V)
    V_forward = V(i)*ones(length(v_star),1);
    vi(:,i) = sqrt(sqrt(V_forward.^4./4+v_star.^4)-V_forward.^2./2); %[m/s] induced velocity in forward flight
end

for i = 1:length(W)
    Pi(i,:) = (W(i).*vi(i,:))./(1000); %[kW] induced power for forward flight
end

%Parassite power
kf = 1.2;
cf = 0.0040;
kw = 2.04;
f = kf*cf*kw*A;
P_p = 1/2*rho_sl*f*V.^3*1e-3; %[kW] parassite power

%Tail rotor power

T_tr = (k*W)./(Omega*x_tr).*sqrt(W/(2*rho_sl*A)); %[N] Tail rotor thrust
DL_tr = T_tr/A_tr; %[N/m^2] disk load for the tail rotor

u0_tr = sqrt(DL_tr/(2*rho_sl)); %[m/s] induced velocity on tail rotor
P_tr = T_tr.*u0_tr*1e-3; %[kW] Tail rotor's (induced) power on forward flight

%Other powers
%By convention we consider the other powers due to: losses due to non
%uniform velocty profile, tip losses and wake downstream the rotor
Pi_star = k*W.*sqrt(W/(2*rho_sl*A))*1e-3; %[kW] induced power (main rotor) in hover
P_oth = 0.17*Pi_star;
for i = 1:length(P_oth)
    P_other(i,:) = P_oth(i)*ones(1,length(V));
end
P_tr = P_tr*ones(1,length(V));
P_tot = Pi+Po+P_p+P_tr+P_other;


hold on
for i = 1:length(W)
    figure(1)
    plot(V,P_tot(i,:),'LineWidth',1.5,'DisplayName',strcat('W  = ',' ',num2str(W(i),5),' N'))

end
title('Requested power for forward flight (@sea level)')
xlabel ('Forward velocity V [m/s]')
ylabel ('Requested Power [kW]')

plot(V,Pa,'LineWidth',1.5,'LineStyle','--','DisplayName','Available power')
legend('show','Location','northwest')

hold off

%la maggior parte della potenza richiesta ad alta velocità è dovuta alla
%potenza indotta

for i = 1:2:length(W)
    hold on
    figure(2)
    plot(V,Pa-P_tot(i,:),'LineWidth',1.5,'DisplayName',strcat('W  = ',' ',num2str(W(i),5),' N'))
end
title('Excess power for forward flight (@sea level)')
xlabel ('Forward velocity V [m/s]')
ylabel ('Excess Power \DeltaP [kW]')
legend('show','Location','best')
hold off

for i = 1:2:length(W)
    hold on
    figure(3)
    plot(V,Pa-Pi(i,:),'LineWidth',1.5,'DisplayName',strcat('W  = ',' ',num2str(W(i),5),' N'))
end
title('Excess induced power for forward flight (@sea level)')
xlabel ('Forward velocity V [m/s]')
ylabel ('Excess Induced Power \DeltaP_i [kW]')
legend('show','Location','best')
hold off

%Po and Pp are invariant with weight
figure(4)
plot(mu,P_tot(1,:),mu,Pi(1,:),mu,Po(1,:),mu,P_p(1,:),mu,P_tr(1,:),mu,P_other(1,:))
legend('tot','ind','prof','paras','tr','oth')

%% Power profile with variable altitude and velocity

W = W_TO; %[N] the operative mass
h = [0:1000:5000]'; %[m] the altitude
T_SL = 288.15; %[K] Temperature at sea level
a = -6.5/1000; %[K/m] temperature gradient
T = T_SL+a.*h; %[K] temperature profile
rho = rho_sl*(T./T_SL).^4.25; %[kg/m^3] air density with altitude
Pa = ne*Pa_sl*bhp_w*eta.*rho./rho_sl.*(1/1000); %[kW] available power for 2 engines


Po = (rho*A*(Omega*R)^3.*Cp0)./(1000); %[kW]profile power with no reverse flow effect

v_star = sqrt(W./A).*sqrt(1./(2*rho)); %[m/s] induced velocity in hover

% for i = 1:length(V)
%     for j = 1:length(rho)
%         vi(j,i) = sqrt(sqrt(V(i).^4./4+v_star(j).^4)-V(i).^2./2) %[m/s] induced velocity in forward flight
%     end
%
% end
for i = 1:length(V)
    V_forward = V(i)*ones(length(v_star),1);
    vi(:,i) = sqrt(sqrt(V_forward.^4./4+v_star.^4)-V_forward.^2./2); %[m/s] induced velocity in forward flight
end

Pi = (W*vi)./(1000); %[kW] induced power for forward flight


%Parassite power
for i = 1:length(rho)
    P_p(i,:) = 1/2*rho(i)*f*V.^3*1e-3; %[kW] parassite power
end
%Tail rotor power

T_tr = (k*W)./(Omega*x_tr).*sqrt(W./(2*rho*A)); %[N] Tail rotor thrust
DL_tr = T_tr/A_tr; %[N/m^2] disk load for the tail rotor

u0_tr = sqrt(DL_tr./(2*rho)); %[m/s] induced velocity on tail rotor
P_tr = T_tr.*u0_tr*1e-3; %[kW] Tail rotor's (induced) power on forward flight

%Other powers
%By convention we consider the other powers due to: losses due to non
%uniform velocty profile, tip losses and wake downstream the rotor
Pi_star = k*W.*sqrt(W./(2*rho*A))*1e-3; %[kW] induced power (main rotor) in hover
P_oth = 0.17*Pi_star;
for i = 1:length(P_oth)
    P_other(i,:) = P_oth(i)*ones(1,length(V));
end

P_tot = Pi+Po+P_p+P_tr+P_other;


for i = 1:length(rho)
    hold on
    figure(4)
    plot(V,P_tot(i,:),'LineWidth',1.5,'DisplayName',strcat('h  = ',' ',num2str(h(i),4),' m'))

end

title('Requested power for forward flight (operative mass)')
xlabel ('Forward velocity V [m/s]')
ylabel ('Requested Power [kW]')
legend('show','Location','northwest')

hold off
%%
%Graphs - Po, Pi, Pp, Ptr, Poth
%Po
for i = 1:length(rho)
    hold on
    figure(5)
    plot(V,Po(i,:),'LineWidth',1.5,'DisplayName',strcat('h  = ',' ',num2str(h(i),4),' m'))
end
title('Profile power for forward flight (operative mass)')
xlabel ('Forward velocity V [m/s]')
ylabel ('P_o [kW]')
legend('show','Location','best')
hold off

%Pi
for i = 1:length(rho)
    hold on
    figure(6)
    plot(V,Pi(i,:),'LineWidth',1.5,'DisplayName',strcat('h  = ',' ',num2str(h(i),4),' m'))
end
title('Induced power for forward flight (operative mass)')
xlabel ('Forward velocity V [m/s]')
ylabel ('P_i [kW]')
legend('show','Location','best')
hold off

%Pp
for i = 1:length(rho)
    hold on
    figure(7)
    plot(V,P_p(i,:),'LineWidth',1.5,'DisplayName',strcat('h  = ',' ',num2str(h(i),4),' m'))
end
title('Parasite power for forward flight (operative mass)')
xlabel ('Forward velocity V [m/s]')
ylabel ('P_p [kW]')
legend('show','Location','best')
hold off

%%
%The deltaP is done taking as the available power that at sea level
%conditions
for i = 1:length(rho)
    dP = Pa(1)*ones(1,length(V))-P_tot(i,:);
    hold on
    figure(8)
    plot(V,dP,'LineWidth',1.5,'DisplayName',strcat('h  = ',' ',num2str(h(i),4),' m'))
end
title('Excess power for forward flight (operative mass)')
xlabel ('Forward velocity V [m/s]')
ylabel ('\DeltaP [kW]')
legend('show','Location','best')
hold off


%UTILE COMMENTARE L'ANDAMENTO DELLE DIVERSE POTENZA CON LA QUOTA E CON IL
%PESO, OLTRE CHE CON LA VELOCITà -> METTERE TUTTO INSIEME -DONE!!!!!

%% Lift-to-Drag ratio

%We still consider W = W operative (takeoff)
% @sl

%Power at sea level (only those useful for the L/D)
P_tot_sl = P_tot(1,:);
Po_sl = Po(1,:);
Pi_sl = Pi(1,:);

L_D  =  (W.*V)./(P_tot_sl*1000); %[-] L/D whole helicopter. N.B.: we multiply the power by 1000 to convert it back to W

L_D_mr = (W.*V)./((Po_sl+Pi_sl)*1000);


figure()
hold on
plot(V,L_D_mr,V,L_D,'LineWidth',2)
title('Lift-to-Drag ratio')
xlabel ('Forward velocity V [m/s]')
ylabel ('L/D')

pos=find(L_D==max(L_D));
V_E_max = V(pos); %[m/s] velocity of main efficiency

plot(V_E_max,max(L_D),'k*')

legend('Main rotor','Helicopter','V_{E_{max}}','Location','best')
hold off

%% Rate of climb R/C
% assuming constant weight = operative weight

% 1 - At sea level
Pa_vec = Pa*ones(1,length(V)); %Available power taken as a vector the sixe of the velocity vector
R_C_sl = ((Pa_vec(1,:)-P_tot_sl)).*1000./W;

figure()
hold on
plot(V,R_C_sl,'LineWidth',2)
title('Rate of climb (@sea level)')
xlabel ('Forward velocity V [m/s]')
ylabel ('R/C [m/s]')

pos = find(R_C_sl==max(R_C_sl));
V_P_min = V(pos); %[m/s] velocity of minimum requested power
R_C_max = R_C_sl(pos); %[m/s] maximum rate of climb

plot(V_P_min,max(R_C_sl),'k*')

hold off

%REMEMBER TO DO R/C at CEILING

%% Ceiling condition
%W = operative weight
clc

h = 0:1:10000; %[m] altitude
V = 0:1:120;
for i=1:length(h)

    for j=1:length(V)

        [Pa, Po, Pi, P_p, P_tr, P_other, P_tot] = power_forwardflight(W, h(i), V(j),1);

        if abs(Pa-P_tot)<70
            dP=abs(Pa-P_tot);
            h_c = h(i); % ceiling altitude
            V_c = V(j); % velocity at ceiling altitude
            
        end
    end

end
h_c
V_c

%to plot the ceiling condition we must evaluate the power at the ceiling
%altitude (now known to us) throughout the velocity profile -> we do the same calculations using  


[Pa, Po, Pi, P_p, P_tr, P_other, P_tot] = power_forwardflight(W,h_c,V);

figure()
plot(V,P_tot,V,Pa,V_c,P_tot(find(V==V_c)),'k*','LineWidth',2)
title(strcat('Power profile @ ceiling condition h_c = ', num2str(h_c), ' m'))
xlabel ('Forward velocity V [m/s]')
ylabel ('Power [kW]')
legend ('Total power requested', 'Available power','Ceiling condition','Location','best')