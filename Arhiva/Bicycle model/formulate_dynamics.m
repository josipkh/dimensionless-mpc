% formulate system dynamics with the bicycle model
A1 = -(2*Cf + 2*Cr)/m/V;
A2 = (2*Cr*Lr - 2*Cf*Lf)/m/V;
A3 = (2*Cr*Lr - 2*Cf*Lf)/Iz/V;
A4 = -(2*Cf*Lf^2 + 2*Cr*Lr^2)/Iz/V;

% the error w.r.t. road model (Rajamani, p.36), second input (road preview) is ignored?
A = [0   1     0    0;
     0   A1  -A1*V  A2;
     0   0     0    1;
     0   A3  -A3*V  A4];
 
% Rajamani, p.30, eq.(2.31)
% A = [0   1    0    0;
%      0   A1   0   A2-V;
%      0   0    0    1;
%      0   A3   0   A4];

B = [0; 2*Cf/m; 0; 2*Lf*Cf/Iz];
C = eye(4);

[Ad,Bd] = c2d(A,B,Ts);

