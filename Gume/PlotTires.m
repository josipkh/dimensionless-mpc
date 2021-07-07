%% Plot a full Magic Formula (5.2), simplified MF and a linear model
clear;clc;close all

% Vehicle parameters
load ../ParamsFull
g = CONST.GRAVITY;
m = VEHICLE.MASS;
L = VEHICLE.WHEEL_BASE;
lf = VEHICLE.LF;
lr = VEHICLE.LR;

% Number of points
nPoints = 200;

% Forward velocity
vx = 10;    % [m/s]

%% Pure longitudinal test case
Fz      = ones(nPoints,1).*m*g*lr/2/L;      % vertical load         (N)
kappa	= linspace(-0.3,0.3, nPoints)';     % longitudinal slip 	(-) (-1 = locked wheel)
alpha	= ones(nPoints,1).*0;               % side slip angle    	(radians)
gamma	= ones(nPoints,1).*0;               % inclination angle 	(radians)
phit 	= ones(nPoints,1).*0;               % turnslip            	(1/m)
Vx   	= ones(nPoints,1).*vx;              % forward velocity   	(m/s)

% Create a string with the name of the TIR file
TIRfile = 'MF_205_60R15_V91.tir';

% Select a Use Mode
useMode = 111;

% Wrap all inputs in one matrix
inputs = [Fz kappa alpha gamma phit Vx];

% Store the output from mfeval in a 2D Matrix
output = mfeval(TIRfile, inputs, useMode);

% Extract the variables from output
Fx = output(:,1);

% Find the linear range coefficient
[~,idx0] = min(abs(kappa));
d1 = diff(Fx);
d2 = diff(kappa);
Cfx = d1(idx0) / d2(idx0);  % derivative at the origin

kappaLinearMax = 5/100;
kappaLinearMin = -kappaLinearMax;
[~,idxMax] = min(abs(kappa-kappaLinearMax));
[~,idxMin] = min(abs(kappa-kappaLinearMin));
kappaLinear = kappa(idxMin:idxMax);
Fx_linear = Cfx * kappaLinear;

% Use a simplified Magic formula
Cx = VEHICLE.TIRE_CX;
Dx = VEHICLE.TIRE_DX;
Ex = VEHICLE.TIRE_EX;
Kx = VEHICLE.TIRE_KX;
Bx = Kx / (Cx*Dx);
Bxk = Bx .* kappa;
mux = Dx*sin(Cx*atan(Bxk-Ex*(Bxk-atan(Bxk))));
Fxs = mux .* Fz;

figure
plot(kappa,Fx./Fz,'b')
hold on
plot(kappaLinear,Fx_linear./Fz(idxMin:idxMax),'r')
plot(kappa,Fxs./Fz,'g')
grid on
title('Longitudinal force')
xlabel('Slip [\%]')
ylabel('Fx/Fz [-]')
legend('MFeval','Linear','Simplified')

%% Pure lateral test case
Fz      = ones(nPoints,1).*m*g*lr/2/L;              % vertical load         (N)
kappa	= ones(nPoints,1).*0;                       % longitudinal slip 	(-) (-1 = locked wheel)
alpha	= linspace(-0.1745,0.1745, nPoints)';       % side slip angle    	(radians)
gamma	= ones(nPoints,1).*0;                       % inclination angle 	(radians)
phit 	= ones(nPoints,1).*0;                       % turnslip            	(1/m)
Vx   	= ones(nPoints,1).*vx;                      % forward velocity   	(m/s)

% Create a string with the name of the TIR file
TIRfile = 'MF_205_60R15_V91.tir';

% Select a Use Mode
useMode = 111;

% Wrap all inputs in one matrix
inputs = [Fz kappa alpha gamma phit Vx];

% Store the output from mfeval in a 2D Matrix
output = mfeval(TIRfile, inputs, useMode);

% Extract the variables from output
Fy = output(:,2);

% Find the linear range coefficient
[~,idx0] = min(abs(alpha));
d1 = diff(Fy);
d2 = diff(alpha);
Cfy = - d1(idx0) / d2(idx0);  % derivative at the origin

alphaLinearMax = deg2rad(4);
alphaLinearMin = -alphaLinearMax;
[~,idxMax] = min(abs(alpha-alphaLinearMax));
[~,idxMin] = min(abs(alpha-alphaLinearMin));
alphaLinear = alpha(idxMin:idxMax);
Fy_linear = Cfy * alphaLinear;

% Use a simplified Magic formula
Cy = VEHICLE.TIRE_CY;
Dy = VEHICLE.TIRE_DY;
Ey = VEHICLE.TIRE_EY;
Ky = VEHICLE.TIRE_KY;
By = Ky / (Cy*Dy);
Bya = By .* alpha;
muy = Dy*sin(Cy*atan(Bya-Ey*(Bya-atan(Bya))));
Fys = muy .* Fz;

figure
plot(rad2deg(alpha),Fy./Fz,'b')
hold on
plot(rad2deg(alphaLinear),-Fy_linear./Fz(idxMin:idxMax),'r')
plot(rad2deg(alpha),Fys./Fz,'g')
grid on
title('Lateral force')
xlabel('Slip angle [deg]')
ylabel('Fy/Fz [-]')
legend('MFeval','Linear','Simplified')