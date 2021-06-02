close all;clc

vx0 = 10;
xi0 = [0; vx0; 0; 0; 250/79*vx0*ones(4,1)];

u0 = 50*ones(4,1);
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