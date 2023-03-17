function [Pa, Po, Pi, Pp, P_tr, P_other, P_tot] = power_forwardflight(W,h,V,c)
%POWER_FORWARDFLIGHT - calculates the requested power [kW] for forward flight considering a forward velocity from 0 to 120 m/s
%   W = weight to use in calculations
%   h = altitude
%   V = forward velocity
%   if 3 input arguments -> result in kW, else result in W

%----------- Input (change desired data for each helicopter) -----------%
g = 9.81; %[m/s^2] gravity @sl
rho_sl = 1.225; %[kg/m^3] density @sl
nb = 4; %blade number
c = 0.391; %[m] mean aerodynamic chord
R = 6.4; %[m] blade radius
Omega = 35.63; %[rad/s] angular velocity - main rotor
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
T_SL = 288.15; %[K] Temperature at sea level
a = -6.5/1000; %[K/m] temperature gradient


%Parassite power coefficients
kf = 1.2;
cf = 0.0040;
kw = 2.04;
f = kf*cf*kw*A;


%Temperature and density
T = T_SL+a.*h; %[K] temperature profile
rho = rho_sl*(T./T_SL).^4.25; %[kg/m^3] air density with altitude


Pa = ne*Pa_sl*bhp_w*eta.*rho./rho_sl*ones(1,length(V)); %[kW] available power for 2 engines


mu = V./(Omega*R); %rapporto di avanzamento
%We shall value the required power without considering the reverse flow
%effect
Cp0 = (sigma*Cd0)/8*(1+K.*mu.^2); %profile drag coefficient without reverse flow (mu<0.2)
Po = (rho*A*(Omega*R)^3.*Cp0); %[kW]profile power with no reverse flow effect

v_star = [sqrt(W./A).*sqrt(1./(2*rho))]; %[m/s] induced velocity in hover

for i = 1:length(V)
    V_forward = V(i)*ones(length(v_star),1);
    vi(:,i) = sqrt(sqrt(V_forward.^4./4+v_star.^4)-V_forward.^2./2); %[m/s] induced velocity in forward flight
end

for i = 1:length(W)
    Pi(i,:) = (W(i).*vi(i,:)); %[kW] induced power for forward flight
end

Pp = 1/2*rho*f.*V.^3; %[kW] parassite power

T_tr = (k*W)./(Omega*x_tr).*sqrt(W/(2*rho*A)); %[N] Tail rotor thrust
DL_tr = T_tr/A_tr; %[N/m^2] disk load for the tail rotor

u0_tr = sqrt(DL_tr/(2*rho)); %[m/s] induced velocity on tail rotor
P_tr = T_tr.*u0_tr*ones(1,length(V)); %[kW] Tail rotor's (induced) power on forward flight

%Other powers
%By convention we consider the other powers due to: losses due to non
%uniform velocty profile, tip losses and wake downstream the rotor
Pi_star = k*W.*sqrt(W/(2*rho*A)); %[kW] induced power (main rotor) in hover
P_oth = 0.17*Pi_star;
for i = 1:length(P_oth)
    P_other(i,:) = P_oth(i)*ones(1,length(V));
end
if nargin==3 %sends result in kW
    Pa = Pa*1e-3;
    Pi = Pi*1e-3;
    Po = Po*1e-3;
    Pp = Pp*1e-3;
    P_tr = P_tr*1e-3;
    P_other = P_other*1e-3;
    P_tot = (Pi+Po+Pp+P_tr+P_other);
   
else %result in W
    P_tot = (Pi+Po+Pp+P_tr+P_other);
    
end


end