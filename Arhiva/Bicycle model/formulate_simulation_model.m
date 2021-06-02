clear;clc;

syms y dy theta dtheta real
syms df real
syms Ts real

%% full size vehicle
vehicle_params = 2;
f = 1;
load_vehicle_parameters;

A1 = -(Cf+Cr)/m;
A2 = (Cr*Lr-Cf*Lf)/m;
A3 = (Cr*Lr-Cf*Lf)/Iz;
A4 = -(Cf*Lf^2+Cr*Lr^2)/Iz;

A = [0   1    0   0;
     0 A1/V -A1 A2/V;
     0   0    0   1;
     0 A3/V -A3 A4/V];
B = [0; Cf/m; 0; Lf*Cf/Iz];

f = A * [y;dy;theta;dtheta] + B * df;

matlabFunction([f;dtheta],'file','OdeBicycleModelLarge','vars',{Ts,[y;dy;theta;dtheta;df]});

%% lab scale vehicle
vehicle_params = 1;
f = 1;
load_vehicle_parameters;

A1 = -(Cf+Cr)/m;
A2 = (Cr*Lr-Cf*Lf)/m;
A3 = (Cr*Lr-Cf*Lf)/Iz;
A4 = -(Cf*Lf^2+Cr*Lr^2)/Iz;

A = [0   1    0   0;
     0 A1/V -A1 A2/V;
     0   0    0   1;
     0 A3/V -A3 A4/V];
B = [0; Cf/m; 0; Lf*Cf/Iz];

f = A * [y;dy;theta;dtheta] + B * df;

matlabFunction([f;dtheta],'file','OdeBicycleModelSmall','vars',{Ts,[y;dy;theta;dtheta;df]});