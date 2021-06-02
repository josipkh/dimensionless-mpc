%% Fit a simplified Magic (5.2) formula to data from a .tir file
clear;clc;close all

% Vehicle parameters
g = 9.81;
m = 1194;
L = 2.59;
lf = 1.01;
lr = 1.58;

% Number of points
nPoints = 200;

% Forward velocity
vx = 10;    % [m/s]

%% Pure longitudinal force
Fz      = ones(nPoints,1).*m*g*lr/2/L;      % vertical load         (N)
kappa	= linspace(-0.3,0.3, nPoints)';        % longitudinal slip 	(-) (-1 = locked wheel)
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

% plot
figure
plot(kappa,Fx./Fz,'b')
hold on
grid on
title('mux-kappa')
xlabel('Slip [\%]')
ylabel('Fx/Fz [-]')

%% Longitudinal fit
F = Fx./Fz;
sx = kappa;

D = 1.2069; Dx = D;
Cs = [1.35:0.02:1.5 1.5];
Es = -4:0.5:-0.5;
Ks = 30:1:100;

K_best = nan;
C_best = nan;
E_best = nan;
error_best = inf;

j=0;
for K = Ks
    for C = Cs
        B = K/(C*D);
        for E = Es
            F_test = D*sin(C*atan(B*sx-E*(B*sx-atan(B*sx))));
            error_test = sum(abs(F-F_test));
            if error_test < error_best
                K_best = K;
                C_best = C;
                E_best = E;
                error_best = error_test;
            end
        end
    end
    j = j+1;
end

Kx = K_best;
Cx = C_best;
Ex = E_best;
Bx = Kx/(Cx*Dx);
B_best = K_best/(C_best*Dx);

Fx_fit = D*sin(C_best*atan(B_best*sx-E_best*(B_best*sx-atan(B_best*sx))));
plot(kappa,Fx_fit,'r')
legend('MFeval','fitted')

%% Pure lateral force
kappa	= ones(nPoints,1).*0;                       % longitudinal slip 	(-) (-1 = locked wheel)
alpha	= linspace(-0.1745,0.1745, nPoints)';       % side slip angle    	(radians)
inputs = [Fz kappa alpha gamma phit Vx];
output = mfeval(TIRfile, inputs, useMode);
Fy = output(:,2);

figure
plot(rad2deg(alpha),Fy./Fz,'b')
hold on
grid on
title('muy-alpha')
xlabel('Slip angle [deg]')
ylabel('Fy/Fz [-]')

%% Lateral fit
F = Fy./Fz;

Dy = 1.1; D = Dy;
Cs = -2:0.1:2;
Es = -5:1:5;
Ks = -100:10:100;

D_best = nan;
C_best = nan;
E_best = nan;
error_best = inf;

j=0;
for K = Ks
    for C = Cs
        B = K/(C*D);
        for E = Es
            F_test = D*sin(C*atan(B*alpha-E*(B*alpha-atan(B*alpha))));
            error_test = sum(abs(F-F_test));
            if error_test < error_best
                K_best = K;
                C_best = C;
                E_best = E;
                error_best = error_test;
            end
        end
    end
    j = j+1;
end

Ky = K_best;
Cy = C_best;
Ey = E_best;
By = Ky/(Cy*Dy);

Fy_fit = D*sin(Cy*atan(By*alpha-Ey*(By*alpha-atan(By*alpha))));
plot(rad2deg(alpha),Fy_fit,'r')
legend('MFeval','fitted')
