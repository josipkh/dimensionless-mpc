p1 = Lf/(Lf+Lr);
p2 = Lr/(Lf+Lr);
p3 = Cf*(Lf+Lr)/m/V^2;
p4 = Cr*(Lf+Lr)/m/V^2;
p5 = Iz/m/(Lf+Lr)^2;

M = diag([Lf+Lr V 1 V/(Lf+Lr)]);
Minv = M\eye(size(A,1));

Ap = (Lf+Lr)/V * (Minv*A*M);
Bp = (Lf+Lr)/V * (Minv*B);
Cp = C*M;

% Ap = [0 1 0 0;
%       0 -p3-p4 p3+p4 p2*p4-p1*p3;
%       0 0 0 1;
%       0 (p2*p4-p1*p3)/p5 (p1*p3-p2*p4)/p5 -(p1^2*p3+p2^2*p4)/p5];
% Bp = [0 p3 0 p1*p3/p5]';

%Ap==(Lf+Lr)/V*M\A*M  % (check difference, 10e-15)
%Bp==(Lf+Lr)/V*M\*B

Tsp = Ts*V/(Lf+Lr);  % normalize time units
[Apd, Bpd] = c2d(Ap,Bp,Tsp);