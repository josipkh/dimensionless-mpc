%% Formulate a simulation model (load transfer, aerodynamics)
clear;clc;

% create symbolic variables for states and inputs
syms acc_long_ms2 acc_lat_ms2 vx vy yaw_rate real
syms deltaf real
syms w_FL w_FR w_RL w_RR real  % wheel speeds
syms Mp_FL Mp_FR Mp_RL Mp_RR real  % propulsion moment on each wheel
syms Ts real  % sampling time, for differential equations (ode45 compatibility)

x = [vx vy yaw_rate w_FL w_FR w_RL w_RR]';  % state vector
u = [Mp_FL Mp_FR Mp_RL Mp_RR]';  % input vector

% define vehicle parameters
load ParamsFull
% chassis
m       = VEHICLE.MASS;
Jz      = VEHICLE.INERTIA_Z;
L       = VEHICLE.WHEEL_BASE;
lf      = VEHICLE.LF;
lr      = VEHICLE.LR;
track_F = VEHICLE.TRACK_FRONT;
track_R = VEHICLE.TRACK_REAR;
hcg     = VEHICLE.COG_HEIGHT;
% actuators
ks      = VEHICLE.STEERING_RATIO;
% wheels and tires
rw      = VEHICLE.WHEEL_RADIUS;
Jw      = VEHICLE.WHEEL_INERTIA;
Cfx     = VEHICLE.SLIP_STIFF;
Crx     = Cfx;             
Cfy     = VEHICLE.CORNERING_STIFF;
Cry     = Cfy;
% aerodynamics
hcp     = hcg;                  % aerodynamic center of pressure
Cd      = VEHICLE.DRAG_COEFF;
A       = VEHICLE.FRONT_AREA;
rho     = CONST.AIR_DENSITY;

g   = CONST.GRAVITY;

% MF_205_60R15_V9
% Magic formula parameters (dry asphalt), fitted to the .tir data
% longitudinal force
Dx = VEHICLE.TIRE_DX;
Cx = VEHICLE.TIRE_CX;
Kx = VEHICLE.TIRE_KX;
Bx = VEHICLE.TIRE_BX;
Ex = VEHICLE.TIRE_EX;
% lateral force
Dy = VEHICLE.TIRE_DY;
Cy = VEHICLE.TIRE_CY;
Ky = VEHICLE.TIRE_KY;
By = VEHICLE.TIRE_BY;
Ey = VEHICLE.TIRE_EY;

disp('Parameters loaded')

%syms Fx_FL Fx_FR Fx_RL Fx_RR real  % from the tire model
%syms Fy_FL Fy_FR Fy_RL Fy_RR real

steer_angle_FL = deltaf / ks;
steer_angle_FR = deltaf / ks;
steer_angle_RL = 0;
steer_angle_RR = 0;

%% aerodynamics

Fx_aero = (vx.^2).*(Cd.*A.*rho)./2;
Fz_aero_f = 0;  %(long_speed_ms.^2).*(Cl_f.*A.*ro)./2;
Fz_aero_r = 0;  %(long_speed_ms.^2).*(Cl_r.*A.*ro)./2;

%% normal forces

Fz_F = (m.*g.*lr + Fz_aero_f.*(lf+lr) - Fx_aero.*hcp - acc_long_ms2.*m.*hcg) ./ (lf + lr) ;
Fz_R = (m.*g.*lf + Fz_aero_r.*(lf+lr) + Fx_aero.*hcp + acc_long_ms2.*m.*hcg) ./ (lf + lr) ;

Fz_FL =  (Fz_F ./ 2) - (acc_lat_ms2.*m.*hcg ./ track_F).*(lr ./ (lf + lr) );
Fz_FR =  (Fz_F ./ 2) + (acc_lat_ms2.*m.*hcg ./ track_F).*(lr ./ (lf + lr) );
Fz_RL =  (Fz_R ./ 2) - (acc_lat_ms2.*m.*hcg ./ track_R).*(lf ./ (lf + lr) );
Fz_RR =  (Fz_R ./ 2) + (acc_lat_ms2.*m.*hcg ./ track_R).*(lf ./ (lf + lr) );

%% tyre translational velocities

% coordinate system axis in same direction as car coordinate system
vx_FL = vx - yaw_rate*track_F/2;
vx_FR = vx + yaw_rate*track_F/2;
vx_RL = vx - yaw_rate*track_R/2;
vx_RR = vx + yaw_rate*track_R/2;

vy_FL = vy + yaw_rate*lf;
vy_FR = vy + yaw_rate*lf;
vy_RL = vy - yaw_rate*lr;
vy_RR = vy - yaw_rate*lr;

% rotated coordinate system through steer angle
vx_w_FL = vx_FL*cos(steer_angle_FL) + vy_FL*sin(steer_angle_FL);
vx_w_FR = vx_FR*cos(steer_angle_FR) + vy_FR*sin(steer_angle_FR);
vx_w_RL = vx_RL*cos(steer_angle_RL) + vy_RL*sin(steer_angle_RL);
vx_w_RR = vx_RR*cos(steer_angle_RR) + vy_RR*sin(steer_angle_RR);

vy_w_FL = -vx_FL*sin(steer_angle_FL) + vy_FL*cos(steer_angle_FL);
vy_w_FR = -vx_FR*sin(steer_angle_FR) + vy_FR*cos(steer_angle_FR);
vy_w_RL = -vx_RL*sin(steer_angle_RL) + vy_RL*cos(steer_angle_RL);
vy_w_RR = -vx_RR*sin(steer_angle_RR) + vy_RR*cos(steer_angle_RR);

%% slip angle

slip_angle_FL = vy_w_FL/vx_w_FL;
slip_angle_FR = vy_w_FR/vx_w_FR;
slip_angle_RL = vy_w_RL/vx_w_RL;
slip_angle_RR = vy_w_RR/vx_w_RR;

%% slip ratio

slip_ratio_FL = (w_FL*rw - vx_w_FL) / vx_w_FL;
slip_ratio_FR = (w_FR*rw - vx_w_FR) / vx_w_FR;
slip_ratio_RL = (w_RL*rw - vx_w_RL) / vx_w_RL;
slip_ratio_RR = (w_RR*rw - vx_w_RR) / vx_w_RR;

%% longitudinal tire forces (simplified Magic formula)

mux_FL = Dx*sin(Cx*atan(Bx*slip_ratio_FL-Ex*(Bx*slip_ratio_FL-atan(Bx*slip_ratio_FL))));
mux_FR = Dx*sin(Cx*atan(Bx*slip_ratio_FR-Ex*(Bx*slip_ratio_FR-atan(Bx*slip_ratio_FR))));
mux_RL = Dx*sin(Cx*atan(Bx*slip_ratio_RL-Ex*(Bx*slip_ratio_RL-atan(Bx*slip_ratio_RL))));
mux_RR = Dx*sin(Cx*atan(Bx*slip_ratio_RR-Ex*(Bx*slip_ratio_RR-atan(Bx*slip_ratio_RR))));
Fx_FL = Fz_FL * mux_FL;
Fx_FR = Fz_FR * mux_FR;
Fx_RL = Fz_RL * mux_RL;
Fx_RR = Fz_RR * mux_RR;

%% lateral tire forces (simplified Magic formula)

muy_FL = Dy*sin(Cy*atan(By*slip_angle_FL-Ey*(By*slip_angle_FL-atan(By*slip_angle_FL))));
muy_FR = Dy*sin(Cy*atan(By*slip_angle_FR-Ey*(By*slip_angle_FR-atan(By*slip_angle_FR))));
muy_RL = Dy*sin(Cy*atan(By*slip_angle_RL-Ey*(By*slip_angle_RL-atan(By*slip_angle_RL))));
muy_RR = Dy*sin(Cy*atan(By*slip_angle_RR-Ey*(By*slip_angle_RR-atan(By*slip_angle_RR))));
Fy_FL = Fz_FL * muy_FL;
Fy_FR = Fz_FR * muy_FR;
Fy_RL = Fz_RL * muy_RL;
Fy_RR = Fz_RR * muy_RR;

%% wheel dynamics

w_dot_FL = (1/Jw)*(Mp_FL  - Fx_FL*rw);
w_dot_FR = (1/Jw)*(Mp_FR  - Fx_FR*rw);
w_dot_RL = (1/Jw)*(Mp_RL  - Fx_RL*rw);
w_dot_RR = (1/Jw)*(Mp_RR  - Fx_RR*rw);

%% state space

ax = (1/m)*(Fx_FL*cos(steer_angle_FL) + Fx_FR*cos(steer_angle_FR) + Fx_RL*cos(steer_angle_RL) + Fx_RR*cos(steer_angle_RR) ...
          - Fy_FL*sin(steer_angle_FL) - Fy_FR*sin(steer_angle_FR) - Fy_RL*sin(steer_angle_RL) - Fy_RR*sin(steer_angle_RR) - Fx_aero) + yaw_rate*vy;

ay = (1/m)*(Fx_FL*sin(steer_angle_FL) + Fx_FR*sin(steer_angle_FR) + Fx_RL*sin(steer_angle_RL) + Fx_RR*sin(steer_angle_RR) ...
          + Fy_FL*cos(steer_angle_FL) + Fy_FR*cos(steer_angle_FR) + Fy_RL*cos(steer_angle_RL) + Fy_RR*cos(steer_angle_RR)) - yaw_rate*vx;

yaw_acc = (1/Jz)*(lf*(Fx_FL*sin(steer_angle_FL) + Fx_FR*sin(steer_angle_FR) + Fy_FL*cos(steer_angle_FL) + Fy_FR*cos(steer_angle_FR)) ...
                - lr*(Fx_RL*sin(steer_angle_RL) + Fx_RR*sin(steer_angle_RR) + Fy_RL*cos(steer_angle_RL) + Fy_RR*cos(steer_angle_RR)) ... 
                + (track_F/2)*(-Fx_FL*cos(steer_angle_FL) + Fx_FR*cos(steer_angle_FR) + Fy_FL*sin(steer_angle_FL) - Fy_FR*sin(steer_angle_FR)) ...
                + (track_R/2)*(-Fx_RL*cos(steer_angle_RL) + Fx_RR*cos(steer_angle_RR) + Fy_RL*sin(steer_angle_RL) - Fy_RR*sin(steer_angle_RR)));

%% Save dynamics equations as MATLAB functions
disp('Formulating state dynamics and simplifying')
f = [ax;ay;yaw_acc;w_dot_FL;w_dot_FR;w_dot_RL;w_dot_RL];
simplify(f);

disp('Writing state dynamics to a function')
matlabFunction([f; u; deltaf; acc_long_ms2; acc_lat_ms2],...
               'file','SimulationModel',...
               'vars',{Ts,[x;u;deltaf;acc_long_ms2;acc_lat_ms2]});
           
disp('Done')           
           
%% Find equilibria
% sol=solve(subs(f,{w_FL,w_FR,w_RL,w_RR,steer_angle_FL,yaw_rate,vy,vx,acc_lat_ms2,acc_long_ms2},...
%                  {31.6456,31.6456,31.6456,31.6456,0,0,0,10,0,0})==0,...
%                  [Mp_FL,Mp_FR,Mp_RL,Mp_RR])
%              
% sol=solve(subs(f(1:4),{steer_angle_FL,yaw_rate,vy,vx,acc_lat_ms2,acc_long_ms2},...
%                  {0,0,0,10,0,0})==0,...
%                  [Mp_FL,Mp_FR,Mp_RL,Mp_RR])
