%Helicopter Flight Mechanics
%Esercitazione 3 - Lynx
%----------- Created by Matheus Padilha -----------%
clear all
close all
clc

%----------- Input -----------%
global Omega gamma mu lambda theta0 theta1s theta1c

Omega = 35.63; %[rad/s] angular velocity - main rotor
gamma = 7.12; %Lock number
mu = 0.2; %advance ratio 
lambda = 0.04; %inflow
theta0 = 0.2; %collective pitch
theta1s = -0.1; %longitudinal cyclic pitch
theta1c = 0.1; %lateral cyclic pitch

p = 20*pi/180; % roll - change as requested
q = 30*pi/180; % pitch - change as requested

disp('Helicopter Flight Mechanics')
disp('Esercitazione 3 - Lynx')
disp(' ')
disp('----Rotor dynamics - Flapping mode----')
disp(' ')

%% Harmonic balance method

beta0 = gamma/8*((1+mu^2)*theta0-4/3*lambda+4/3*mu*theta1s+2/3*mu*p/Omega);
beta1s = -1/(1+mu^2/2)*(4/3*mu*beta0-q/Omega-16/gamma*p/Omega)+theta1c;
beta1c = (-2*mu*(4/3*theta0-lambda)-(1+3/2*mu^2)*theta1s+(16/gamma*q/Omega-p/Omega))/(1-mu^2/2);

beta_harmonic_deg = [beta0, beta1s beta1c]*180/pi;
disp('Harmonic Balance Method')
disp('beta0 [deg]')
disp(rad2deg(beta0))
disp('beta1s [deg]')
disp(rad2deg(beta1s))
disp('beta1c [deg]')
disp(rad2deg(beta1c))

%% Galerkin's Method
nstep = 1000; %number of integration steps
dpsi = (2*pi)/nstep;
psi = 0;
A = zeros(5,5);
B = zeros(5,1);
for i=1:nstep
    psi = psi + dpsi;
    [a,b] = galerkin_method_flap(psi,p,q);
    A = A+a;
    B = B+b;
end

beta = A\B;
beta_deg = beta'*180/pi