%=== Saab_9-3_data ========================================================
% Numerical vehicle data representative for a typical Saab 9-3 Sedan
%	Ref: 2003 Saab 9-3 Production No. 1787
%   Created by Matthijs Klomp 2006-08-08
%   Not modified, it's a real vehicle!
%==========================================================================

g=9.81;             %[m/s^2]  Gravitational constant

%% --- Measured Performance --- http://www.saab.com -----------------------
% http://www.saab.com/main/GLOBAL/en/model/93_S/techspecs.shtml

topspeed=220;             %[km/h]    Top speed (unconstrained)
acc_0_100=8.5;            %[s] Acceleration  0 to 100 km/h
acc_60_100=8.5;           %[s] Acceleration 60 to 100 km/h (4:th gear)
acc_80_120=12.5;          %[s] Acceleration 80 to 120 km/h (5:th gear)
fuel_cons=[10.8 5.7 7.6]; %[l/100km] Fuel consumption City/Highway/Combined 
%                                    according to 1999/100 EC directive.

%% --- Chassis ------------------------------------------------------------
%
% Configuration: Front wheel drive, transverse engine
%
m=1675;             %[kg]     Curb weight (full tank, no driver or pass.)
m_unspr=185;        %[kg]     Unsprung mass all four wheels.
m_maxload=450;      %[kg]     Maximum load (driver not included)
Ixx = 540;          %[kg*m^2] Total mass moment of inertia, around CoG
Iyy = 2398;         %[kg*m^2] Total mass moment of inertia, around CoG
Izz = 2617;         %[kg*m^2] Total mass moment of inertia, around CoG
I(1,1)=Ixx; I(2,2)=Iyy; I(3,3)=Izz; % Inertia matrix
L=2.675;            %[m]      Wheel base
l_f=0.4*L;          %[m]      Distance along X-axis from CoG to front axle, MODIFY WEIGHT HERE!!!!
l_r=L-l_f;          %[m]      Distance along X-axis from CoG to rear axle
B_f=1.517;          %[m]      Track width front
B_r=1.505;          %[m]      Track width rear
h=0.543;            %[m]      Distance along Z-axis from CoG to ground.
A_f=2.17;           %[m^2]    Front area
c_d=0.30;           %[-]      Air drag coefficient 
C_zf=30.8e3;        %[N/m]    Front Suspension Ride Rate (w/o tire)/ side**
C_zr=29.9e3;        %[N/m]    Rear Suspension Ride rate (w/o tire)/ side**
D_zf=4500;          %[Ns/m]   Suspension damping front, per side**
D_zr=3500;          %[Ns/m]   Suspension damping rear, per side**
C_af=16.22;         %[N/rad]  Front auxiliary roll stiffness (w/o springs)
C_ar=7.837;         %[N/rad]  Rear  auxiliary roll stiffness (w/o springs)
hr_f=0.045;         %[m]      Front roll center height above ground.
hr_r=0.101;         %[m]      Rear  roll center height above ground.
                    %         ** If bicycle model is used, multiply by 2
                 
%% --- Hook and Pulley ----------------------------------------------------
l_h=1.115;          %[m]      Distance from rear axle to hook for pulley
m2=1500;            %[kg]     Pulley total mass
Iz2=500;            %[kg*m^2] Pulley mass moment of inertia, around CoG *
l_p=2.5;            %[m]      Distance from hook to axle on the pulley *
l_m=2.0;            %[m]      Distance from hook to CoG on the pulley *

%% --- Tire / Wheel -------------------------------------------------------
% Type: Michelin Primacy 215/55 R16
C_tire = [-0.00235  30.7]; % Tire stiffness parameters
% To calculate the cornering stiffness: Calpha = C_tire(1)*Fz^2 + C_tire(2)*Fz
Calphaf=9.3e4;      %[N/rad]  Cornering stiffness / tire @ front axle load
Calphar=7.5e4;      %[N/rad]  Cornering stiffness / tire @ rear  axle load
Calphap=5.0e4;      %[N/rad]  Cornering stiffness for one tire on pulley
deltaT=20;          %[Nm]     Rolling resistance torque for one axle
RRconst=0.0164;     %[-]      Rolling resistance constant
RRvsqd=4.0e-8;	    %[-]      Rolling resistance velocity squared term
C_zt=220e3;	        %[N/m]    Tire stiffnes for one tire.
Ty_const=1.89;	    %[m]      Tire lag to lateral force builup
Tx_const=0.12;      %[m]      Tire lag to longitudianl force builup
radius_wheel=0.316; %[m]      Radius of wheel
J_wheel=2;          %[kg*m^2] Moment of inertia for two wheels and axle

%% --- Environment * ------------------------------------------------------
air=1.3;            %[kg/m^3] Air density

%% --- Driveline ----------------------------------------------------------
%
% Engine configuration: Inline 4 cyl, 16 valves DOHC, Turbo
% Displacement: 1998 cm^3 
% Max power:  129 kW (175 hp) @ 5500 rpm
% Max torque: 265 Nm 	      @ 2500 - 4000 rpm
% Gearbox: Manual 5 speed

%% Ratio for each gear (1-5) [-] ------- Manual 5-speed F35 PR2
ratio_vec=[3.385 1.76 1.179 0.894 0.66];
%ratio_vec=3.5:-0.001:0.5; %%CVT gearbox
ratio_rev=-3.166;   %[-]      Ratio of reverse gear
ratio_final=4.05;   %[-]      Ratio of final gear

% Engine speed vector corresponding to torque vector
w_vec_rpm = [ 0 500 750 1000 1250 1500 1750 2000 2500 ...
             3000 3500 4000 4500 5000 5500 6000 6500];  % [rpm]
w_vec     = w_vec_rpm*pi/30;                            % [rad/s]  

% Engine maximum torque vector (throttle = 100%) [Nm]
Tmax_vec=[  0.0000  72.3550 108.6000 144.8500 175.8800 206.8800 ...
          234.8900 251.9300 264.9800 265.0000 265.0000 263.0000 ...
          255.0200 240.0300 223.0300 204.0400 187.0300];

% Engine minimum torque vector (throttle = 0%) [Nm]
Tmin_vec=[0 0 -5 -12 -20 -28 -36.5 -44.5 -56.8 ...
             -64 -65 -65 -65 -65 -65 -65 -65];

%% --- Steering -----------------------------------------------------------
ratio_steer = 15.9;    %[-]      Ratio of steering gear

%% --- Brake system -------------------------------------------------------
% Brake moment gain per unit master cylinder pressure / wheel
BrakeGainFront = 293;  %[Nm/Mpa] 
BrakeGainRear  = 111;  %[Nm/Mpa] 

%% --- End of file --------------------------------------------------------
