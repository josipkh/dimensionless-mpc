clear;clc;close all

% Vehicle parameters
g = 9.81;
m = 1194;
L = 2.59;
lf = 1.01;
lr = 1.58;
vx = 10;    % [m/s]

% Pure longitudinal test case - front tire
Fz      = m*g*lr/2/L;      % vertical load         (N)
kappa	= 0.2;             % longitudinal slip 	(-) (-1 = locked wheel)
alpha	= 0;               % side slip angle    	(radians)
gamma	= 0;               % inclination angle 	(radians)
phit 	= 0;               % turnslip            	(1/m)
Vx   	= vx;              % forward velocity   	(m/s)

% Create a string with the name of the TIR file
TIRfile = 'MF_205_60R15_V91.tir';

% Select a Use Mode
useMode = 111;

tic
for n = 1:10
% Wrap all inputs in one matrix
inputs = [Fz kappa alpha gamma phit Vx];

% Store the output from mfeval in a 2D Matrix
output = mfeval(TIRfile, inputs, useMode);
Fx = output(1);
n
end
toc