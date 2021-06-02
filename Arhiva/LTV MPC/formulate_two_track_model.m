%% Formulate two-track model dynamics
clear;clc

% define vehicle parameters
% full size vehicle parameters from Saab_93_datasheet.m
m=1675;         %[kg]     Curb weight (full tank, no driver or pass.)
Jz=2617;        %[kg*m^2] Total mass moment of inertia, around CoG
l=2.675;        %[m]      Wheel base
lf=0.4*l;       %[m]      Distance along X-axis from CoG to front axle
lr=l-lf;        %[m]      Distance along X-axis from CoG to rear axle
w=1.517/2;      %[m]      Track width front/2
rw=0.316;       %[m]      Radius of wheel
Jw=1;           %[kg*m^2] Moment of inertia for one wheel and half axle
Cfx=3e4;        %[N]      Tire stiffness parameters
Crx=3e4;        %(Fx=Cx*sx, from the compendium)
Cfy=9.3e4;      %[N/rad]  From the datasheet
Cry=9.3e4;

syms m lf lr w Jz Jw rw real  % vehicle and wheel parameters
syms Cfx Cfy Crx Cry real  % tire coefficients

% create symbolic variables
syms deltaf real  % wheel steering angle, external input
syms Ts real  % sampling time, for difference equations
syms vy vx theta thetad wfl wfr wrl wrr deltaf real  % states
syms Tfl Tfr Trl Trr real  % inputs
xi = [vy vx theta thetad wfl wfr wrl wrr]';
u = [Tfl Tfr Trl Trr]';

% wheel velocities (vehicle frame)
vflx = vx - w * thetad;  vfly = vy + lf * thetad;
vfrx = vx + w * thetad;  vfry = vy + lf * thetad;
vrlx = vx - w * thetad;  vrly = vy - lr * thetad;
vrrx = vx + w * thetad;  vrry = vy - lr * thetad;

% wheel velocities (wheel frame)
vwflx = vflx * cos(deltaf) + vfly * sin(deltaf);   vwrlx = vrlx;
vwfrx = vfrx * cos(deltaf) + vfry * sin(deltaf);   vwrrx = vrrx;
vwfly = -vflx * sin(deltaf) + vfly * cos(deltaf);  vwrly = vrly;
vwfry = -vfrx * sin(deltaf) + vfry * cos(deltaf);  vwrry = vrry;

% tire slips and slip angles
sflx = rw * wfl / vwflx - 1;    alphafl = atan2(vwfly, vwflx);
sfrx = rw * wfr / vwfrx - 1;    alphafr = atan2(vwfry, vwfrx);
srlx = rw * wrl / vwrlx - 1;    alpharl = atan2(vwrly, vwrlx);
srrx = rw * wrr / vwrrx - 1;    alpharr = atan2(vwrry, vwrrx);

% linear tire force model
Fwflx = Cfx * sflx;  Fwfly = -Cfy * alphafl;
Fwfrx = Cfx * sfrx;  Fwfry = -Cfy * alphafr;
Fwrlx = Crx * srlx;  Fwrly = -Cry * alpharl;
Fwrrx = Crx * srrx;  Fwrry = -Cry * alpharr;

% longitudinal and lateral forces
Fflx = Fwflx * cos(deltaf) - Fwfly * sin(deltaf);   Frlx = Fwrlx;
Ffrx = Fwfrx * cos(deltaf) - Fwfry * sin(deltaf);   Frrx = Fwrrx;
Ffly = Fwflx * sin(deltaf) + Fwfly * cos(deltaf);   Frly = Fwrly;
Ffry = Fwfrx * sin(deltaf) + Fwfry * cos(deltaf);   Frry = Fwrry;

% vehicle dynamics
x = thetad * vy + 1 / m * (Fflx + Ffrx + Frlx + Frrx);
y = -thetad * vx + 1 / m * (Ffly + Ffry + Frly + Frry);
z = 1 / Jz * (lf * (Ffly + Ffry) - lr * (Frly + Frry) + w * (-Fflx + Ffrx - Frlx + Frrx));
w1 = 1 / Jw * (Tfl - rw * Fwflx);
w2 = 1 / Jw * (Tfr - rw * Fwfrx);
w3 = 1 / Jw * (Trl - rw * Fwrlx);
w4 = 1 / Jw * (Trr - rw * Fwrrx);

f = [y; x; thetad; z; w1; w2; w3; w4];
f = simplify(f);

% Jacobians
DfDx = jacobian(f,xi);
DfDu = jacobian(f,u);
%% Save dynamics equations as MATLAB functions
matlabFunction(f,'file','twoTrackModel','vars',{xi,u,deltaf});
matlabFunction(DfDx,DfDu,'file','twoTrackModelJacobian','vars',{xi,u,deltaf});
matlabFunction([f;u;deltaf],'file','odeTwoTrackModel','vars',{Ts,[xi;u;deltaf]});
save fDxDu.mat f DfDx DfDu
%% Save vehicle parameters
params.m=m;
params.Jz=Jz;
params.l=l;
params.lf=lf;
params.lr=lr;
params.w=w;  
params.rw=rw;
params.Jw=Jw;
params.Cfx=Cfx;
params.Cfy=Cfy;
params.Crx=Crx;
params.Cry=Cry;
save params.mat params
%% Find equilibria
%sol=solve(subs(f,{deltaf,thetad,Trl,Trr,Tfl,Tfr,vy,vx},{0,0,0,0,0,0,0,10})==0,[wfl,wfr,wrl,wrr])
%sol=solve(subs(f,{deltaf,thetad,Trl,Trr,Tfl,Tfr,vy},{0,0,0,0,0,0,0})==0,[vx,wfl,wfr,wrl,wrr])
%sol=solve(subs(f,{deltaf,thetad,Trl,Trr},{0,0,0,0})==0,[xi([1:2,4:end]);u])