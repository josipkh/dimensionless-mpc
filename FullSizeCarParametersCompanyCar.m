%==============================================================================================
% Numerical vehicle data representative for the Company Car vehicle from IPG Carmaker examples
%==============================================================================================
clear;clc;

%% --- Chassis ------------------------------------------------------------
VEHICLE.MASS                = 1600;         %[kg]       Total mass
VEHICLE.INERTIA_Z           = 2393.665;     %[kg*m^2]   Total mass moment of inertia, around CoG
VEHICLE.WHEEL_BASE          = 2.622;        %[m]        Distance between the axles
VEHICLE.LF                  = 1.311;        %[m]        Distance along X-axis from CoG to front axle
VEHICLE.LR                  = 1.311;        %[m]        Distance along X-axis from CoG to rear axle
VEHICLE.TRACK_FRONT         = 2.622;        %[m]        Track width front
VEHICLE.TRACK_REAR          = 2.622;        %[m]        Track width rear
VEHICLE.COG_HEIGHT          = 0.526;        %[m]        Distance along Z-axis from CoG to ground
VEHICLE.FRONT_AREA          = 2.156;        %[m^2]      Front area
VEHICLE.DRAG_COEFF          = 0.37;         %[-]        Air drag coefficient

%% --- Tire / Wheel -------------------------------------------------------
% Type front and rear: MF_205_60R15_V91
% estimated from the .tir file, Fx=Cx*sx, Fy=-Cy*alphay
VEHICLE.SLIP_STIFF          = 116377;       %[N]        Longitudinal slip stiffness
VEHICLE.CORNERING_STIFF     = 72705;        %[N/rad]    Cornering stiffness
VEHICLE.WHEEL_RADIUS        = 0.318;        %[m]        Radius of wheel
VEHICLE.WHEEL_INERTIA       = 1.22;         %[kg*m^2]   Moment of inertia for one wheel and half-axle

% tire linear range:
VEHICLE.LIN_TIRE_KAPPA_MIN  = -0.3;         %[-]        Linear range, minimal slip
VEHICLE.LIN_TIRE_KAPPA_MAX  =  0.3;         %[-]        Linear range, maximal slip
VEHICLE.LIN_TIRE_ALPHA_MIN  = -0.0436;      %[rad]      Linear range, minimal slip angle, 2.5 deg
VEHICLE.LIN_TIRE_ALPHA_MAX  =  0.0436;      %[rad]      Linear range, maximal slip angle, 2.5 deg

% simplified Magic formula parameters:
% Fx = Fz * D*sin(C*atan(B*sx-E*(B*sx-atan(B*sx))))
% longitudinal force
VEHICLE.TIRE_DX             = 1.2069;
VEHICLE.TIRE_CX             = 1.35;
VEHICLE.TIRE_KX             = 30;
VEHICLE.TIRE_BX             = VEHICLE.TIRE_KX/(VEHICLE.TIRE_CX*VEHICLE.TIRE_DX);
VEHICLE.TIRE_EX             = -0.5;
% lateral force
VEHICLE.TIRE_DY             = 1.1;
VEHICLE.TIRE_CY             = -1.2;
VEHICLE.TIRE_KY             = -20;
VEHICLE.TIRE_BY             = VEHICLE.TIRE_KY/(VEHICLE.TIRE_CY*VEHICLE.TIRE_DY);
VEHICLE.TIRE_EY             = 0;

%% --- Environment --------------------------------------------------------
CONST.AIR_DENSITY           = 1.205;        %[kg/m^3]   Air density (default in Carmaker)
CONST.GRAVITY               = 9.81;         %[m/s^2]    Gravitational constant
CONST.GROUND_FRICTION       = 1.0;          %[-]        Friction coefficient between the ground and the tires

%% --- Driveline ----------------------------------------------------------
VEHICLE.MAX_POWER           = 20;           %[kW]       Maximal engine power
VEHICLE.MAX_TORQUE          = 250;          %[Nm]       Maximal engine torque
VEHICLE.MAX_TORQUE_RATE     = 1000;         %[Nm/s]     Maximal engine torque rate
VEHICLE.MAX_ENG_SPEED       = 8000;         %[rpm]      Maximal engine speed

%% --- Steering -----------------------------------------------------------
VEHICLE.STEERING_RATIO      = 13;           %[-]        Steering wheel angle / road wheel angle
VEHICLE.MAX_STEERING_RATE   = 2*pi;         %[rad/s]    Maximal steering angle rate
VEHICLE.MAX_STEERING_ANGLE  = 2*pi;         %[rad]      Maximal steering wheel angle

%% Save the parameters in a .mat file
save ParamsFull VEHICLE CONST
disp('Full size car parameters saved')
