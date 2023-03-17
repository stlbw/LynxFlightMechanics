%Helicopter Flight Mechanics
%Esercitazione 2 - Lynx
%----------- Created by Matheus Padilha -----------%
clear
close all
clc
%----------- Input -----------%
g = 9.81; %[m/s^2] gravity @sl
rho_sl = 1.225; %[kg/m^3] density @sl
W = 4313.7*g; %[N] operative mass
Cd0  =  1/100; %profile drag coefficient
a = 2*pi; %CL_alpha
k = 1/90; %induced power factor

nb = 4; %blade number
c = 0.391; %[m] mean aerodynamic chord
R = 6.4; %[m] blade radius
A = pi*R^2;
sigma = (nb*c*R)/A;%solidity ratio
Omega = 35.63; %[rad/s] angular velocity - main rotor

R_tr = 1.106; %[m] tail rotor radius
A_tr = pi*R_tr^2; %[m^2] tail rotor's area
x_tr = 7.660; %[m] distance from main rotor to tail rotor
nb_tr = 4; %blade number - tail rotor
Cd0_tr = Cd0;
sigma_tr = 0.208;%solidity ratio tail rotor - table
Omega_tr = 6*Omega; %[rad/s] tail rotor's angular velocity
%c_tr = R_tr/7; %[m] mean aerodynamic chord - tail rotor -> decided not to
%use it as we have directly the solidility ratio

hR = 1.274; %[m] height of main rotor hub above fuselage reference point
hR_tr = 1.146; %[m] height of tail rotor hub above fuselage reference point
ltR = 7.66; %[m] distance of tail rotor hub aft of fuselage reference point
lR = 0; %[m] distance of main rotor hub aft of fuselage reference point
fR = 0; %[m] distance of main rotor hub to the fuselage reference point on lateral axis direction

disp('Helicopter Flight Mechanics')
disp('Esercitazione 2 - Lynx')
disp(' ')
disp('----Equilibrium states - Hovering @sea level----')
disp(' ')

%% Longitudinal plane
h = 100; %[m] height above sea level
rho = rho_sl;
V=0; %[m/s] forward speed
alpha = 0; % AoA, hover condition

% inflow parameter
lambda_i = 1/(Omega*R).*sqrt(W./(2.*rho*A)); %inflow parameter in hovering
C_T = lambda_i.^2*2; %thrust coefficient

lambda = lambda_i - V./(Omega*R).*sin(alpha); %inflow parameter for forward flight -> in this case it will be the same as lambda_i as we are indeed in hovering

% collective pitch theta0
% from the thrust equation from the Blade Element theory, assuming T = W: 

theta0 = 3/2*(lambda + 4.*W./(nb.*rho*a*c*Omega^2*R^3)); %[rad] collective pitch angle 

Theta = 0; % pitch angle - trim condition -> demonstrate it in the report!

a_1s = 0; % (B1 = a1 -> pilot maintains the cyclic command in a fixed position)

%% Latero-directional plane
Cd = Cd0 + k*a^2*(theta0 - lambda)^2;
Q = 1/8*nb*(Cd+lambda*(a*(theta0-lambda))).*rho*c*Omega^2*R^4; %[N.m] the reaction torque due to main rotor -> Blade Element Theory

T_tr = Q./ltR; %[N] the thrust needed from the tail rotor to balance the reaction torque from the main rotor



b_1s = -(T_tr*hR_tr)/(W*hR); %[rad] angle between control axis NFP and hub plane HB

Phi = -T_tr/W - b_1s; %[rad] roll angle 

lambda_t = sqrt(T_tr/(2*rho*A_tr*Omega_tr^2*R_tr^2)); %inflow parameter for the tail rotor
C_T_tr = lambda_t^2*2;

c_tr = (sigma_tr*A_tr)/(nb_tr*R_tr); %[m] tail rotor blade's aerodynamic mean chord

theta0_tr =  3/2*(lambda_t + 4.*T_tr./(nb_tr.*rho*a*c_tr*Omega_tr^2*R_tr^3)); %[rad] tail rotor's collective pitch angle
%% Display parameters
disp('Longitudinal plane')
disp(' ')
disp(' Inflow parameter [-]:')
disp(lambda)
disp(' Collective pitch Theta0- main rotor [deg]:')
disp(rad2deg(theta0))
disp(' Pitch angle Theta [deg]:')
disp(rad2deg(Theta))
disp(' a_1s [deg]:')
disp(rad2deg(a_1s))
disp(' ')
disp('------------------------------------------------')
disp(' ')
disp('Latero-directional plane')
disp(' ')
disp(' Angle between NFP and HB b_1s[deg]:')
disp(rad2deg(b_1s))
disp(' Roll angle Phi [deg]:')
disp(rad2deg(Phi))
disp(' Inflow parameter - tail rotor [-]:')
disp(lambda_t)
disp(' Collective pitch - tail rotor [deg]:')
disp(rad2deg(theta0_tr))