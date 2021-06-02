%% Bicycle model from (2000)
yalmip('clear')
clear

% parameters
m = 1.47;  % kg
Iz = 0.024;  % kgm^2
Lf = 0.13;  % m
Lr = 0.15;  % m
Cf = 1.53;  % N/deg
Cr = 1.53;  % N/deg
V = 3;  % m/s

A1 = -(Cf+Cr)/m;
A2 = (Cr*Lr-Cf*Lf)/m;
B1 = Cf/m;
B2 = Cr/m;
A3 = (Cr*Lr-Cf*Lf)/Iz;
A4 = -(Cf*Lf^2+Cr*Lr^2)/Iz;
B3 = Lf*Cf/Iz;
B4 = -Lr*Cr/Iz;

% system dynamics
A = [0   1    0   0;
     0 A1/V -A1 A2/V;
     0   0    0   1;
     0 A3/V -A3 A4/V];
B = [0   0;
     B1 B2;
     0   0;
     B3 B4];
 
nx = 4;  % no. of states
nu = 2;  % no. of inputs

x0 = [0;0];  % initial condition