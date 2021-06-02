clear;clc;close all

% load simulation data
s1 = load('Simulation results/acc_pudu_osqp/acc_pudu_osqp.mat');
s2 = load('Simulation results/acc_udu_osqp/acc_udu_osqp.mat');

% compare the states
figure
plot((s1.Xis-s2.Xis)')
legend({'vy [m/s]','vx [m/s]','\theta [rad]','d\theta [rad/s]',...
           '\omega FL [rad/s]','\omega FR [rad/s]','\omega RL [rad/s]','\omega RR [rad/s]'})
title('State difference')

% compare the inputs
figure
stairs((s1.Us-s2.Us)')
legend('T_{fl}','T_{fr}','T_{rl}','T_{rr}')
title('Input difference [Nm]')

% compare the input rates
figure
plot((diff(s1.Us,1,2)-diff(s2.Us,1,2))')
legend('\Delta T_{fl}','\Delta T_{fr}','\Delta T_{rl}','\Delta T_{rr}')
title('Input rate difference [Nm]')
