%% Formulate two-track model dynamics (prediction model)
% TODO: load transfer, Ackermann & steering ratio, aerodynamics
clear;clc

% define vehicle parameters
% full size vehicle parameters for examples/demo_compact in Carmaker
g = 9.81;               %[m/s^2]  Gravitational constant
% chassis
m=1194;                 %[kg]     Total mass
Jz=1645.687;            %[kg*m^2] Total mass moment of inertia, around CoG
l=2.59;                 %[m]      Wheel base
lf=1.01;                %[m]      Distance along X-axis from CoG to front axle
lr=1.58;                %[m]      Distance along X-axis from CoG to rear axle
w=1.522/2;              %[m]      Track width front/2
% actuators
Tmax=500;               %[Nm]     Input torque upper bound
Tmin=-Tmax;             %[Nm]     Input torque lower bound
dTmax=1000;             %[Nm/s]   Input torque rate upper bound
dTmin=-dTmax;           %[Nm/s]   Input torque rate lower bound
% wheels and tires
rw=0.318;               %[m]      Radius of wheel
Jw=1.22;                %[kg*m^2] Moment of inertia for one wheel and half axle
Cfx=103050;             %[N]      Tire stiffness parameters, estimated from the .tir file
Crx=103050;             %(Fx=Cx*sx, from the compendium)
Cfy=67043;              %[N/rad]  MF_205_60R15_V91
Cry=67043;
% simplified Magic formula parameters
Dx = 1.2069;
Cx = 1.35;
Kx = 30;
Bx = Kx/(Cx*Dx);
Ex = -0.5;

Dy = 1.1;
Cy = -2;
Ky = -20;
By = Ky/(Cy*Dy);
Ey = 1;

% tire linear range
slipMax=3/100;          %[-]      Allowed slip upper bound (for stability)
slipMin=-slipMax;       %[-]      Allowed slip lower bound (for stability)
alphaMax=deg2rad(3.4);  %[rad]    Allowed slip angle upper bound (for stability)
alphaMin=-alphaMax;     %[rad]    Allowed slip angle lower bound (for stability)

% syms m lf lr w Jz Jw rw real  % vehicle and wheel parameters
% syms Cfx Cfy Crx Cry real  % tire coefficients
% syms slipMax slipMin alphaMax alphaMin real  % for state constraints

% create symbolic variables
syms deltaf real  % wheel steering angle, external input
syms Ts real  % sampling time, for difference equations
syms vy vx thetad wfl wfr wrl wrr deltaf real  % states
syms Tfl Tfr Trl Trr real  % inputs
xi = [vy vx thetad wfl wfr wrl wrr]';
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

% % simplified Magic formula
% Fzf = m*g*lr/l/2;
% Fzr = m*g*lf/l/2;
% Fz = [Fzf; Fzf; Fzr; Fzr];
% sxs = [sflx; sfrx; srlx; srrx];
% alphas = [alphafl; alphafr; alpharl; alpharr];
% 
% muxs = Dx*sin(Cx*atan(Bx*sxs-Ex*(Bx*sxs-atan(Bx*sxs))));
% muys = Dy*sin(Cy*atan(By*alphas-Ey*(By*alphas-atan(By*alphas))));
% 
% Fwflx = muxs(1) * Fz(1);
% Fwfrx = muxs(2) * Fz(2);
% Fwrlx = muxs(3) * Fz(3);
% Fwrrx = muxs(4) * Fz(4);
% 
% Fwfly = muys(1) * Fz(1);
% Fwfry = muys(2) * Fz(2);
% Fwrly = muys(3) * Fz(3);
% Fwrry = muys(4) * Fz(4);

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

disp('Formulating state dynamics and simplifying')
f = [x; y; z; w1; w2; w3; w4];
f = simplify(f);

% Jacobians
disp('Calculating the Jacobians')
DfDx = jacobian(f,xi);
DfDu = jacobian(f,u);

% %% Save dynamics equations as MATLAB functions
% disp('Writing the Jacobians to a function')
% matlabFunction(DfDx,DfDu,'file','TwoTrackModelJacobian','vars',{xi,u,deltaf});
% disp('Writing state dynamics to a function')
% matlabFunction([f;u;deltaf],'file','OdeTwoTrackModel','vars',{Ts,[xi;u;deltaf]});
% disp('Done')
% % matlabFunction(f,'file','twoTrackModel','vars',{xi,u,deltaf});
% % save fDxDu.mat f DfDx DfDu

%% Save vehicle parameters
% params.m=m;
% params.Jz=Jz;
% params.l=l;
% params.lf=lf;
% params.lr=lr;
% params.w=w;  
% params.rw=rw;
% params.Jw=Jw;
% params.Cfx=Cfx;
% params.Cfy=Cfy;
% params.Crx=Crx;
% params.Cry=Cry;
% params.Tmax=Tmax;
% params.Tmin=Tmin;
% params.dTmax=dTmax;
% params.dTmin=dTmin;
% params.slipMax=slipMax;
% params.slipMin=slipMin;
% params.alphaMax=alphaMax;
% params.alphaMin=alphaMin;
% 
% params.nx = length(xi);
% params.nu = length(u);
% 
% save Params.mat params

%% Find equilibria
%sol=solve(subs(f,{deltaf,thetad,Trl,Trr,Tfl,Tfr,vy,vx},{0,0,0,0,0,0,0,10})==0,[wfl,wfr,wrl,wrr])
%sol=solve(subs(f,{deltaf,thetad,Trl,Trr,Tfl,Tfr,vy},{0,0,0,0,0,0,0})==0,[vx,wfl,wfr,wrl,wrr])
%sol=solve(subs(f,{deltaf,thetad,Trl,Trr},{0,0,0,0})==0,[xi([1:2,4:end]);u])
