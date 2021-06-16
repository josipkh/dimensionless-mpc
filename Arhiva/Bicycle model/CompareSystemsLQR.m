%% parameters IRS 2001
m = 6.02;  % kg
Iz = 0.153;  % kgm^2
Lf = 0.137;  % m
Lr = 0.222;  % m
Cf = 43*2;  % N/rad (2 tires)
Cr = 43*2;  % N/rad (2 tires)
V = 1.98;  % m/s
Ts = 0.01;

%% parameters full size (Carmaker demo compact)
m = 1194;  % kg
Iz = 1645.687;  % kgm^2
Lf = 1.01;  % m
Lr = 1.58;  % m
Cf = 67043*2;  % N/rad (2 tires)
Cr = 67043*2;  % N/rad (2 tires)
V = 15;  % m/s
Ts = 0.01;

%% formulate system dynamics with the bicycle model
A1 = -(Cf+Cr)/m;
A2 = (Cr*Lr-Cf*Lf)/m;
A3 = (Cr*Lr-Cf*Lf)/Iz;
A4 = -(Cf*Lf^2+Cr*Lr^2)/Iz;

A = [0   1    0   0;
     0 A1/V -A1 A2/V;
     0   0    0   1;
     0 A3/V -A3 A4/V];
B = [0; Cf/m; 0; Lf*Cf/Iz];
C = eye(4);

[Ad,Bd] = c2d(A,B,Ts);

%% formulate Pi-groups
M = diag([Lf+Lr V 1 V/(Lf+Lr)]);
Minv = M\eye(size(A,1));

Ap = (Lf+Lr)/V * (Minv*A*M);
Bp = (Lf+Lr)/V * (Minv*B);
Cp = C*M;

Tsp = Ts*V/(Lf+Lr);  % normalize time units
[Apd, Bpd] = c2d(Ap,Bp,Tsp);

p1 = Lf/(Lf+Lr);
p2 = Lr/(Lf+Lr);
p3 = Cf*(Lf+Lr)/m/V^2;
p4 = Cr*(Lf+Lr)/m/V^2;
p5 = Iz/m/(Lf+Lr)^2;
p = [p1 p2 p3 p4 p5];

%% design LQRs and compare closed-loop poles
[K,S,CLP] = lqr(A,B,diag([1 0 0 0]),1);
[Kp,Sp,CLPp] = lqr(Ap,Bp,diag([1 0 0 0]),1);