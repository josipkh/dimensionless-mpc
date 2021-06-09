% A script to compare nonlinear models with the linearized versions
close all;clc; clear

Ts = 0.05;

vx0 = 10;  % m/s
xi0 = [vx0; 0; 0; 250/79*vx0*ones(4,1)];

u0 = 50*ones(4,1);
[A,B] = TwoTrackModelJacobian(xi0,u0,0);
[Ad,Bd] = c2d(A,B,Ts);

tend = 0.5;
yl = xi0;
for i=1:tend/0.05
    yl = [yl Ad*yl(:,end)+Bd*u0];
end

[t,y]=ode45(@OdeTwoTrackModel,[0,tend],[xi0;u0;0]);
figure;
plot(0.05*[0:tend/0.05],yl(1,:))
hold on
plot(t,y(:,1))
legend('linearized','nonlinear')