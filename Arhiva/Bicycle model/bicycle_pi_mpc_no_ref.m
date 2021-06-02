%% Bicycle model (2001) & Pi groups - MPC controller
clear;clc;close all

% load vehicle parameters
vehicle_params = 1;  % 1: IRS(2001), 2: full size(2000,p.260)
f = 1;  % scaling factor
load_vehicle_parameters;

% system dynamics
Ts = 0.05;  % sampling time
formulate_dynamics;  % A,B,C,Ad,Bd

% Pi-groups
formulate_pi_groups;  % M,Ap,Bp,Cp,Apd,Bpd
% M = eye(4); Apd = Ad; Bpd = Bd;

%% MPC
yalmip('clear');close all;clc

% MPC setup
nx = size(A,1);  % no. of states
nu = size(B,2);  % no. of inputs

Q = 10^6*diag([1 0 1e-3 0]);  % state weights
Q = M*Q*M;
R = 1e-6;  % input weights
N = 10;  % horizon length

% YALMIP data
u = sdpvar(nu*ones(1,N),ones(1,N));
x = sdpvar(nx*ones(1,N+1),ones(1,N+1));

% constraints and cost
constraints = [];
objective = 0;
for k = 1:N
    objective = objective + x{k}'*Q*x{k} + u{k}'*R*u{k};
    constraints = [constraints, x{k+1} == Apd*x{k} + Bpd*u{k}];
    constraints = [constraints, -0.35 <= u{k}<= 0.35];  % ,-5<=x{k+1}(i)<=5
end

solutions_out = {[u{:}], [x{:}]};
parameters_in = x{1};
x = [0.1;0;0;0];  % initial state

controller = optimizer(constraints, objective, ...
    sdpsettings('solver','osqp'),parameters_in,solutions_out);

clf; hold on
n = 100;  % simulation length
xhist = nan(nx,n+1);
xhist(:,1) = x;
uhist = nan(nu,n);
Xhist = nan(nx,n+1);
Xhist(:,1) = M\x;
for i = 1:n
    x = Minv*x;  % convert init. state to pi-space
    inputs = {x};
    [solutions,diagnostics] = controller{inputs};    
    U = solutions{1};
    X = solutions{2};
    if diagnostics == 1
        error('The problem is infeasible');
    end    
    
    subplot(2,1,1);stairs(i:i+length(U)-1,U,'r');title('Predicted input')
    subplot(2,1,2);cla;stairs(i:i+N,X(1,:),'b');
    hold on;
    stairs(1:i,Xhist(1,1:i),'g');title('Prediction+reference+state')
    
    Xhist(:,i+1) = x;
    x = M*x;  % convert back to "normal" space
    x = x + Ts*(A*x + B*U(1));  % very simple simulation
    
    xhist(:,i+1) = x;
    uhist(:,i) = U(1);
    pause(0.1)   
end

figure; t = (1:n)*Ts;
subplot(3,1,1); plot(t,uhist); title('Implemented input')
subplot(3,1,2); plot([0 t],xhist(1,:),'b'); title('Lateral position')
subplot(3,1,3); plot([0 t],xhist(3,:)); title('Yaw angle')