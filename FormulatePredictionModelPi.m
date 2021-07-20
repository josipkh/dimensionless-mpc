%% Formulate two-track model dynamics (prediction model)
clear;clc;
load ParamsFull

% define vehicle parameters
% chassis
m   = VEHICLE.MASS;
Jz  = VEHICLE.INERTIA_Z;
lf  = VEHICLE.LF;
lr  = VEHICLE.LR;
w   = VEHICLE.TRACK_FRONT;
% actuators
Tmax= VEHICLE.MAX_TORQUE;
Tmin=-VEHICLE.MAX_TORQUE;
ks  = VEHICLE.STEERING_RATIO;
% wheels and tires
rw  = VEHICLE.WHEEL_RADIUS;
Jw  = VEHICLE.WHEEL_INERTIA;
Cfx = VEHICLE.SLIP_STIFF;
Crx = Cfx;             
Cfy = VEHICLE.CORNERING_STIFF;
Cry = Cfy;

g   = CONST.GRAVITY;

disp('Parameters loaded')

% syms m lf lr w Jz Jw rw real  % vehicle and wheel parameters
% syms Cfx Cfy Crx Cry real  % tire coefficients
% syms slipMax slipMin alphaMax alphaMin real  % for state constraints

% create symbolic variables
syms deltaf real  % wheel steering angle, external input
syms Ts real  % sampling time, for difference equations
syms vx vy thetad wfl wfr wrl wrr deltaf real  % states
syms Tfl Tfr Trl Trr real  % inputs
xi = [vx vy thetad wfl wfr wrl wrr]';
u = [Tfl Tfr Trl Trr]';

% steering wheel angle to wheel angles
deltafl = deltaf / ks;
deltafr = deltaf / ks;

% wheel velocities (vehicle frame)
vflx = vx - w/2 * thetad;  vfly = vy + lf * thetad;
vfrx = vx + w/2 * thetad;  vfry = vy + lf * thetad;
vrlx = vx - w/2 * thetad;  vrly = vy - lr * thetad;
vrrx = vx + w/2 * thetad;  vrry = vy - lr * thetad;

% wheel velocities (wheel frame)
vwflx =  vflx * cos(deltafl) + vfly * sin(deltafl);   vwrlx = vrlx;
vwfrx =  vfrx * cos(deltafr) + vfry * sin(deltafr);   vwrrx = vrrx;
vwfly = -vflx * sin(deltafl) + vfly * cos(deltafl);  vwrly = vrly;
vwfry = -vfrx * sin(deltafr) + vfry * cos(deltafr);  vwrry = vrry;

% tire slips and slip angles
sflx = rw * wfl / vwflx - 1;    alphafl = atan2(vwfly, vwflx);
sfrx = rw * wfr / vwfrx - 1;    alphafr = atan2(vwfry, vwfrx);
srlx = rw * wrl / vwrlx - 1;    alpharl = atan2(vwrly, vwrlx);
srrx = rw * wrr / vwrrx - 1;    alpharr = atan2(vwrry, vwrrx);

% simplified Magic formula
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
Fflx = Fwflx * cos(deltafl) - Fwfly * sin(deltafl);   Frlx = Fwrlx;
Ffrx = Fwfrx * cos(deltafr) - Fwfry * sin(deltafr);   Frrx = Fwrrx;
Ffly = Fwflx * sin(deltafl) + Fwfly * cos(deltafl);   Frly = Fwrly;
Ffry = Fwfrx * sin(deltafr) + Fwfry * cos(deltafr);   Frry = Fwrry;

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
%% Rewrite the dynamics in dimensionless form
% define new dimensionless variables
syms vxp vyp thetadp wflp wfrp wrlp wrrp Tflp Tfrp Trlp Trrp real
xip = [vxp vyp thetadp wflp wfrp wrlp wrrp]';
up = [Tflp Tfrp Trlp Trrp]';

% define reference values from the system's parameters
speed_ref = sqrt(Cfx*(lf+lr)/m);
ang_speed_ref = sqrt(Cfx/m/(lf+lr));
torque_ref = (lf+lr)*Cfx;
time_ref = sqrt(m*(lf+lr)/Cfx);

% replace the physical variables on the RHS with (reference * dimensionless variable)
disp('Rewriting in dimensionless form')
fp = subs(f,{vx,vy,thetad,wfl,wfr,wrl,wrr,Tfl,Tfr,Trl,Trr},...
            {speed_ref*vxp,speed_ref*vyp,ang_speed_ref*thetadp,...
             ang_speed_ref*wflp,ang_speed_ref*wfrp,ang_speed_ref*wrlp,ang_speed_ref*wrrp,...
             torque_ref*Tflp,torque_ref*Tfrp,torque_ref*Trlp,torque_ref*Trrp});
         
% multiply the RHS with normalization inverses of the LHS to get dx*/dt*
fp = fp * time_ref;
fp(1:2) = fp(1:2) / speed_ref;
fp(3:end) = fp(3:end) / ang_speed_ref;

%% Calculate the Jacobians
disp('Calculating the Jacobians')
DfpDxip = jacobian(fp,xip);
DfpDup = jacobian(fp,up);

%% Save dynamics equations as MATLAB functions
disp('Writing the Jacobians to a function')
matlabFunction(DfpDxip,DfpDup,'file','PredictionModelJacobianPi','vars',{xip,up,deltaf});
disp('Writing state dynamics to a function')
matlabFunction([fp;up;deltaf],'file','PredictionModelPi','vars',{Ts,[xip;up;deltaf]});
disp('Done')
% matlabFunction(f,'file','twoTrackModel','vars',{xi,u,deltaf});

%% Find equilibria
%sol=solve(subs(f,{deltaf,thetad,Trl,Trr,Tfl,Tfr,vy,vx},{0,0,0,0,0,0,0,10})==0,[wfl,wfr,wrl,wrr])
%sol=solve(subs(f,{deltaf,thetad,Trl,Trr,Tfl,Tfr,vy},{0,0,0,0,0,0,0})==0,[vx,wfl,wfr,wrl,wrr])
%sol=solve(subs(f,{deltaf,thetad,Trl,Trr},{0,0,0,0})==0,[xi([1:2,4:end]);u])
