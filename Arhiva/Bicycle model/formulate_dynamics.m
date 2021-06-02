% formulate system dynamics with the bicycle model

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