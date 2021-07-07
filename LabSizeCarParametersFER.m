%==============================================================================================
% Numerical vehicle data representative of the FER lab vehicle
%==============================================================================================
clear;clc;

%% --- Chassis ------------------------------------------------------------
VEHICLE.MASS                = 1.173 + 1;    %[kg]       Total mass; extra mass at CoG to match Pi groups
VEHICLE.INERTIA_Z           = 0.03373;      %[kg*m^2]   Total mass moment of inertia, around CoG
VEHICLE.WHEEL_BASE          = 0.256;        %[m]        Distance between the axles
VEHICLE.LF                  = 0.115;        %[m]        Distance along X-axis from CoG to front axle
VEHICLE.LR                  = 0.141;        %[m]        Distance along X-axis from CoG to rear axle
VEHICLE.TRACK_FRONT         = 0.163;        %[m]        Track width front
VEHICLE.TRACK_REAR          = 0.163;        %[m]        Track width rear
VEHICLE.COG_HEIGHT          = 0.045;        %[m]        Distance along Z-axis from CoG to ground
VEHICLE.FRONT_AREA          = 0;            %[m^2]      Front area
VEHICLE.DRAG_COEFF          = 0;            %[-]        Air drag coefficient

%% --- Tire / Wheel -------------------------------------------------------
% Type front and rear: MF_205_60R15_V91
% estimated from the .tir file, Fx=Cx*sx, Fy=-Cy*alphay
VEHICLE.SLIP_STIFF          = 25;           %[N]        Longitudinal slip stiffness
VEHICLE.CORNERING_STIFF     = 8.25;         %[N/rad]    Cornering stiffness
VEHICLE.WHEEL_RADIUS        = 0.0325;       %[m]        Radius of wheel
VEHICLE.WHEEL_INERTIA       = 3.697e-5;     %[kg*m^2]   Moment of inertia for one wheel and half-axle

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
VEHICLE.TIRE_CY             = -2;
VEHICLE.TIRE_KY             = -20;
VEHICLE.TIRE_BY             = VEHICLE.TIRE_KY/(VEHICLE.TIRE_CY*VEHICLE.TIRE_DY);
VEHICLE.TIRE_EY             = 1;

%% --- Environment --------------------------------------------------------
CONST.AIR_DENSITY           = 1.205;        %[kg/m^3]   Air density (default in Carmaker)
CONST.GRAVITY               = 9.81;         %[m/s^2]    Gravitational constant
CONST.GROUND_FRICTION       = 1.578;        %[-]        Friction coefficient between the ground and the tires

%% --- Driveline ----------------------------------------------------------
VEHICLE.MAX_POWER           = 20;           %[kW]       Maximal engine power
VEHICLE.MAX_TORQUE          = 0.05;         %[Nm]       Maximal engine torque
VEHICLE.MAX_TORQUE_RATE     = 0.2;          %[Nm/s]     Maximal engine torque rate
VEHICLE.MAX_ENG_SPEED       = 8000;         %[rpm]      Maximal engine speed

%% --- Steering -----------------------------------------------------------
VEHICLE.STEERING_RATIO      = 1;            %[-]        Steering wheel angle / road wheel angle
VEHICLE.MAX_STEERING_RATE   = 0.262;        %[rad/s]    Maximal steering angle rate
VEHICLE.MAX_STEERING_ANGLE  = 0.262;        %[rad]      Maximal steering wheel angle

%% Save the parameters in a .mat file
save ParamsLab VEHICLE CONST
disp('Lab scale car parameters saved')