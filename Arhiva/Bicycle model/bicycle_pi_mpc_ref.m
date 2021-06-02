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

%% MPC
yalmip('clear');close all;clc

nx = size(A,1);  % no. of states
nu = size(B,2);  % no. of inputs
Cref = [1 0 0 0];  % tracked states (only lat. position)
Cref = Cref*M;
ny = size(Cref,1);  % no. of tracked states

Q = 10^6*diag([1 0 1e-3 0]);  % state weights
Q = M*Q*M;
%Q = 10^6*eye(nx);  % match state dimensions
R = 1e-6;  % input weights
N = 20;  % horizon length

% YALMIP data
u = sdpvar(nu*ones(1,N),ones(1,N));
x = sdpvar(nx*ones(1,N+1),ones(1,N+1));
r = sdpvar(ny*ones(1,N+1),ones(1,N+1));
ops = sdpsettings('verbose',2);  % print output

% constraints and cost
constraints = [];
%pastu = sdpvar(1);
%constraints = [-.1 <= diff([pastu u{:}]) <= .1];  % constrain du
objective = 0;
for k = 1:N
    objective = objective + (Cref*x{k}-r{k})'*(Cref*x{k}-r{k}) + u{k}'*u{k};
    constraints = [constraints, x{k+1} == Apd*x{k} + Bpd*u{k}];
    constraints = [constraints, -0.35 <= u{k}<= 0.35];  % ,-5<=x{k+1}(i)<=5
end
objective = objective + (C*x{N+1}-r{N+1})'*(C*x{N+1}-r{N+1});  % add terminal cost?

solutions_out = {[u{:}], [x{:}]};
parameters_in = {x{1},[r{:}]};
x = [0;0;0;0];  % initial state

controller = optimizer(constraints, objective, ...
    sdpsettings('solver','gurobi'),parameters_in,solutions_out);

clf; hold on
n = 100;  % simulation length
xhist = nan(nx,n+1);
xhist(:,1) = x;
uhist = nan(nu,n);
rhist = nan(ny,n);
Rhist = nan(ny,n);  % reference in pi-space
Xhist = nan(nx,n+1);  % state in pi-space
Xhist(:,1) = M\x;

% sine reference parameters
Aref = 0.2;
w = 2;
Arefp = Aref/M(1);  % dimensionless
wp = w*(Lf+Lr)/V;  % dimensionless
for i = 1:n
    x = M\x;  % convert init. state to pi-space
    future_r = Aref*sin(w*(i:i+N)*Ts); 
    % convert the reference to pi-space (needs to be generalized)
    future_R = Arefp*sin(wp*(i:i+N)*Tsp);
    
    inputs = {x,future_R};
    [solutions,diagnostics] = controller{inputs};    
    U = solutions{1};
    X = solutions{2};
    if diagnostics == 1
        error('The problem is infeasible');
    end    
    
    subplot(2,1,1);stairs(i:i+length(U)-1,U,'r');title('Predicted input')
    subplot(2,1,2);cla;stairs(i:i+N,X(1,:),'b');
    hold on;stairs(i:i+N,future_R(1,:),'k');
    stairs(1:i,xhist(1,1:i)/M(1),'g');title('Prediction+reference+state')
    
    Xhist(:,i+1) = x;
    Rhist(:,i) = future_R(:,1);
    
    x = M*x;  % convert back to "normal" space
    x = x + Ts*(A*x + B*U(1));
    
    xhist(:,i+1) = x;
    uhist(:,i) = U(1);
    rhist(:,i) = future_r(:,1);
    pause(0.01)   
end

figure; t = (1:n)*Ts; tp = (1:n)*Tsp;
subplot(3,1,1); plot(t,uhist); title('Implemented input')
subplot(3,1,2); plot([0 t],xhist(1,:),'b'); title('Lateral position')
hold on; plot(t,rhist(1,:),'r--');
subplot(3,1,3); plot([0 t],xhist(3,:)); title('Yaw angle')
figure;plot(tp,Rhist,'r--');
hold on;plot([0 tp],Xhist(1,:),'b');title('Pi-space')
