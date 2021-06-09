% A script for comparing dimensional and dimensionless dynamics models
clear;clc;close all

% vehicle parameters
load ParamsFull.mat
m   = VEHICLE.MASS;
L   = VEHICLE.WHEEL_BASE;
Cfx = VEHICLE.SLIP_STIFF;

nx = 7;
nu = 4;

% Pi-groups
Mxi = diag([sqrt(Cfx*L/m),...
           sqrt(Cfx*L/m),...
           sqrt(Cfx/m/L)*ones(1,5)]);
Mxinv = Mxi\eye(nx);
Mu = L*Cfx * eye(nu);
Muinv = Mu\eye(nu);

% approx. steady state as xi0
vx0 = 10;
xi0 = [vx0; 0; 0; 250/79*vx0*ones(4,1)];

% control input
u0 = VEHICLE.MAX_TORQUE / 5 *ones(4,1);

% simulate the system
tend = 5;

[t,y]=ode45(@OdeTwoTrackModel,[0,tend],[xi0;u0;0]);
plot(t,y(:,1))

hold on
[tp,yp]=ode45(@OdeTwoTrackModelPi,[0,tend],[Mxinv*xi0;Muinv*u0;0]);
plot(t,(Mxi(1)*y(:,1)))
legend

%%

[A,B] = twoTrackModelJacobian(xi0,u0,0);
[Ad,Bd] = c2d(A,B,Ts);

tend = 0.5;
yl = xi0;
for i=1:tend/0.05
    yl = [yl Ad*yl(:,end)+Bd*u0];
end

[t,y]=ode45(@odeTwoTrackModel,[0,tend],[xi0;u0;0]);
figure;
plot(0.05*[0:tend/0.05],yl(2,:))
hold on
plot(t,y(:,2))
legend('linearized','nonlinear')

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
