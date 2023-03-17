%Helicopter Flight Mechanics
%Esercitazione 5A e 5B - Lynx
%----------- Created by Matheus Padilha -----------%
clear
close all
clc
%----------- Input -----------%
g = 9.81; % [m/s^2]
delta_theta0 = 1e-2; %[m] collective pitch command
B1 = 1e-2; %[m] longitudinal cyclic pitch command


disp('Helicopter Flight Mechanics')
disp('Esercitazione 5A - Lynx')
disp(' ')
disp('----Helicopter Dynamic response on longitudinal plane (time domain)----')
disp(' ')

%% Uncoupled planes - dynamic response on longitudinal plane - Numeric simulation

[A,~,~, B] = stability_derivatives_Lynx(0);

[eigvec, eigval] = eig(A);

C = eye(4);
D = zeros(4,2);

sys = ss(A,B,C,D);
[P,Z] = pzmap(sys); %returns the poles and the zeros of the system

figure()
hold on
pzmap(sys) %plots the poles of the system (eigenvalues)
grid on

sys1 = tf(sys); %transfer function

%PHUGOID MODE
wn_ph = sqrt(real(eigval(2,2))^2 + imag(eigval(2,2))^2); %natural frequency [rad/s]
zeta = -real(eigval(2,2))/wn_ph; %damping - phugoid mode (zeta<0 == unstable mode)

%% Uncoupled planes - dynamic response on longitudinal plane - Approximate model for phugoid mode
% for the approximate model we consider q_dot = 0 and theta_dot = q
% x_approx = [u, q, theta]'
% u = [theta, B1]'

A_phug = [A(1,1) 0 -g
    A(3,1) A(3,3) 0
    0 1 0];

B_phug = [B(1,1) B(1,2)
    B(3,1) B(3,2)
    0 0];

C_phug = eye(3);
D_phug = zeros(3,2); % 2 == command input

sys_phug = ss(A_phug,B_phug,C_phug,D_phug);

[P_phug, Z_phug] = pzmap(sys_phug); %poles and zeros of the approximated phugoid model


pzmap(sys_phug)
grid on
hold off
legend('Sistema longitudinale', 'Trattazione semplificata')
[eigvec_phug_approx, eigval_phug_approx] = eig(A_phug);

wn_ph_approx = sqrt(real(eigval_phug_approx(2,2))^2+imag(eigval_phug_approx(2,2))^2);
zeta_ph_approx = -real(eigval_phug_approx(2,2))/wn_ph_approx;

disp('------------------------------------------------------------------')
disp('                         PHUGOID MODE (PH)                        ')
disp('Natural frequency [rad/s]')
disp(wn_ph)
disp('Damping')
disp(zeta)
disp('------------------------------------------------------------------')
disp('')
disp('------------------------------------------------------------------')
disp('             PHUGOID MODE - APPROXIMATE SOLUTION (PH)             ')
disp('Natural frequency [rad/s]')
disp(wn_ph_approx)
disp('Damping')
disp(zeta_ph_approx)
disp('------------------------------------------------------------------')

%% Response to the collective pitch command
% NUMERICAL SOLUTION
%We already have the transfer function for the longitudinal plane system.

t = 0:0.1:100; % simulation's temporal evolution
u = delta_theta0*ones(size(t)); %the input command -> collective pitch == 1 cm
w_num = lsim(sys1(2,1),u,t);%N.B.: sys1(2,1) means from input 1 (theta0) to output 2 (w)
figure()
hold on
plot(t,w_num,'b','LineWidth',1.5)

% Stability Augmentation System
% In order to make the system stable, we must have every mode (on
% longitudinal plane in this case) stable. In our case we shall stabilize
% the phugoid mode as it is the only unstable mode, while making sure we do
% not destabilize the other modes.

%The stability augmentation system is added to the longitudinal cyclic
%pitch B1 = B1c + B1_SAS;

K_aug = [0 0 0 0
        0 0 0 0.1];

A_aug = A+B*K_aug;

eig(A_aug);
A_aug = A + B*K_aug;
B_aug = B;
C_aug = C;
D_aug = D;

sys_aug = ss(A_aug,B_aug,C_aug,D_aug);
sys_aug1 = tf(sys_aug);

w_num_aug = lsim(sys_aug1(2,1),u,t);%N.B.: sys1(2,1) means from input 1 (theta0) to output 2 (w)

plot(t,w_num_aug,'c','LineWidth',1.5)

% ANALYTIC SOLUTION
Zw = A(2,2);
Z_theta0 = B(2,1);

w_analytic = -Z_theta0/Zw*delta_theta0*(1-exp(Zw*t));

plot(t,w_analytic,'r','LineWidth',1.5)
hold off
legend('Numerical solution','Augmented system', 'Analytical solution')
title('w(t) dynamic response to \Delta\theta_0')
ylabel('w(t) [m/s]')
xlabel('time [s]')

% Response to the cyclic pitch command
% we keep the temporal evolution as before
t2 = 0:0.1:20;

%NUMERICAL
u2 = B1*ones(size(t2)); %the input command
q_num = lsim(sys1(3,2),u2,t2);
figure()
hold on
plot(t2,q_num,'b','LineWidth',1.5)


q_num_aug = lsim(sys_aug1(3,2),u2,t2);%N.B.: sys1(2,1) means from input 1 (theta0) to output 2 (w)

plot(t2,q_num_aug,'c','LineWidth',1.5)


%ANALYTICAL
M_B1 = B(3,2);
Mq = A(3,3);

q_analytic = -M_B1/Mq*B1*(1-exp(Mq*t2));
plot(t2,q_analytic,'r','LineWidth',1.5)
legend('Numerical solution','Augmented system', 'Analytical solution','Location','best')
title('q(t) dynamic response to \Delta B_1 (\Delta \theta_{1s})')
ylabel('q(t) [rad/s]')
xlabel('time [s]')
hold off


%% plot esercitazione 

long = lsim(sys1(:,1), u, t);
figure()
plot(t,long(:,1))
xlabel('time [s]')
ylabel('u [m/s]')
title('Risposta dinamica del sistema al comando \Delta\theta_0 = 1 cm')
figure()
plot(t,long(:,2))
xlabel('time [s]')
ylabel('w [m/s]')
title('Risposta dinamica del sistema al comando \Delta\theta_0 = 1 cm')
figure()
plot(t,long(:,3))
xlabel('time [s]')
ylabel('q [rad/s]')
title('Risposta dinamica del sistema al comando \Delta\theta_0 = 1 cm')
figure()
plot(t,long(:,4))
xlabel('time [s]')
ylabel('\theta [rad/s]')
title('Risposta dinamica del sistema al comando \Delta\theta_0 = 1 cm')



%% ANALISI RISPOSTA IN FREQUENZA - 5B
%Taken from Prof. Elisa Capello 

disp('Esercitazione 5B - Lynx')
disp(' ')
disp('----Helicopter Dynamic response on longitudinal plane (frequency domain)----')
disp(' ')

[num,den] = tfdata(sys_aug1(4,2),'v');
[z,p,k] = zpkdata(sys_aug1(4,2),'v');

% Poli del sistema aumentato (stabilizzato)
p_aug = eig(A_aug);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Aggiunta di un attuatore
% Attuatore
wn   = 52;      % [rad/s]
zita = 0.7;

% Dinamica attutatore in termini di FDT
act = tf(wn^2,[1 2*zita*wn wn^2]);

% State-space formulation
[A_act,B_act,C_act,D_act] = ssdata(act);

% Matrice aumentata
A_c = [A_aug zeros(4,2);
    zeros(2,4) A_act];

% Autovalori del sistema stabilizzato più attuatore
p_act = eig(A_c);
wn_ph2   =  sqrt(real(p_act(2))^2+imag(p_act(2))^2);
zeta_ph2 = -real(p_act(2))/wn_ph2;

disp('------------------------------------------------------------------')
disp(' PHUGOID MODE (PH) - Augmented and corrected with actuator system ')
disp('Natural frequency [rad/s]')
disp(wn_ph2)
disp('Damping')
disp(zeta_ph2)
disp('------------------------------------------------------------------')

%% Bode di theta
u = 1e-2*ones(size(t2));
% Risposta theta numerica
teta_numerica = lsim(sys1(4,2),u,t2);
% Risposta theta stabilizzata

figure()
grid on
plot(t2,teta_numerica,'b','LineWidth',1.5)
hold on
% Stabilizzazione del sistema
hold on
teta_aug=lsim(sys_aug1(4,2),u,t2); %teta stabilizzata
plot(t2,teta_aug,'c','linewidth',1.5)
hold on
grid on
ylabel('Amplitude')
xlabel('Time [s]')
title('Risposta di \theta al gradino unitario di \Delta B1')
legend('Soluzione numerica','soluzione stabilizzata')

% Aggiunta attuatore
sys_act_teta = sys_aug1(4,2)*act;

figure()
bode(-sys_aug1(4,2),'g')
hold on
bode(-sys_act_teta,'b')
grid on
title('Diagramma di Bode di \Theta/B1')
legend('Sistema stabilizzato', 'Sistema stabilizzato con attuatore', 'Location','best')

%% Normativa ADS-33
%Short period response
omega_180 = 8.52; %[rad/s]
delta_phi = 201; %[deg]
tau_p = deg2rad(delta_phi)/(2*omega_180);

%Mid-term response
%For the mid-term response we consider phugoid's frequency and damping in
%the AUGMENTED AND CORRECTED SYSTEM (with actuator)
x_value = -zeta_ph2*wn_ph2
y_value = wn_ph2*sqrt(1-zeta_ph2^2)


%% Risposta nel tempo r-theta_tr
%%%% SOLO PER CAPIRE RISPOSTA NEL TEMPO PIANO LATERO-DIREZIONALE
Nttr = -0.008;
Nr   = -0.362;
Delta_tr = 1;

a = Nttr/Nr*Delta_tr;
b = Nr;
t2 = 0:0.1:40;
u = ones(size(t2));

r_analitica = a*(1-exp(b*t2));

figure()
plot(t2,r_analitica,'r');
xlabel('Time [s]')
ylabel('r [rad/s]')
grid on

