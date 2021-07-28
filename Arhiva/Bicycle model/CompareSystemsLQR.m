clear;clc;
f = 1;
%% lab size vehicle
vehicle_params = 1;
load_vehicle_parameters;
Ts = 0.01;

formulate_dynamics;
formulate_pi_groups;
p_lab = [p1 p2 p3 p4 p5]

[Klab,Slab,CLPlab] = lqr(A,B,diag([1 0 0 0]),1);
[Kplab,Splab,CLPplab] = lqr(Ap,Bp,diag([1 0 0 0]),1);

% figure; step(ss(Ap,Bp,C,0))
eig(Ap)
%% full size vehicle
vehicle_params = 2;
load_vehicle_parameters;
Ts = 0.01;

formulate_dynamics;
formulate_pi_groups;
p_full = [p1 p2 p3 p4 p5]

[Kfull,Sfull,CLPfull] = lqr(A,B,diag([1 0 0 0]),1);
[Kpfull,Spfull,CLPpfull] = lqr(Ap,Bp,diag([1 0 0 0]),1);

% figure; step(ss(Ap,Bp,C,0))
eig(Ap)