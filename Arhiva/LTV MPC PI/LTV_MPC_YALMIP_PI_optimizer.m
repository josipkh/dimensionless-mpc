%% Two-track model LTV MPC
clear;yalmip('clear');close all;clc
load params.mat  % vehicle parameters
nx = 8;  % no. of states
nu = 4;  % no. of inputs
C = [eye(4),zeros(4)];
ny = size(C,1);  % no. of tracked states

Q = diag([10 100 1e4 2e4]);  % tracked state weights
QN = Q;  % terminal cost weights
R = 1e-3*eye(nu);  % input weights
p = 100;  % slack variable weights

Ts = 0.05;  % sampling time
N = 10;  % prediction horizon
M = 3;  % control horizon (blocking)

% YALMIP data
Ad = sdpvar(nx,nx,'full');
Bd = sdpvar(nx,nu);
xihat = sdpvar(nx*ones(1,N+1),ones(1,N+1));
u = sdpvar(nu*ones(1,M),ones(1,M));  % inputs
xi = sdpvar(nx*ones(1,N+1),ones(1,N+1));  % states
xiref = sdpvar(ny*ones(1,N+1),ones(1,N+1));  % reference
e = sdpvar(16*ones(1,N+1),ones(1,N+1));  % soft constraint variables
pastu = sdpvar(4,1);  % for input rate limits
deltaf = sdpvar(1,1);  % steering angle
ops = sdpsettings('verbose',1,...
                  'solver','osqp',...
                  'osqp.eps_abs',1e-6,...
                  'osqp.eps_rel',1e-6,...
                  'debug',1,...
                  'cachesolvers',1);  % print output

% QP formulation
slipMax = 3.8/100;
slipAngleMax = deg2rad(3.4);
constraints = [];
constraints = [constraints, -500 <= [u{:}] <= 500];  % input constraint
constraints = [constraints, -50 <= diff([pastu u{:}]) <= 50];  % input rate constraint
constraints = [constraints, [e{:}] >= 0];  % slack variables are positive
objective = 0;
for k = 1:N
    objective = objective + (C*xi{k}-xiref{k})'*Q*(C*xi{k}-xiref{k});
    if k > M
        objective = objective + u{M}'*R*u{M};
        constraints = [constraints, xi{k+1} == Ad*xi{k} + Bd*u{M} + xihat{k+1} - Ad*xihat{k} - Bd*pastu];
    else
        objective = objective + u{k}'*R*u{k};
        constraints = [constraints, xi{k+1} == Ad*xi{k} + Bd*u{k} + xihat{k+1} - Ad*xihat{k} - Bd*pastu];
    end
    % state constraints
    % tire slip constraints
    slipConstraints = [    
        params.rw*xi{k+1}(5) - (slipMax + 1)*(cos(deltaf)*(xi{k+1}(2) - xi{k+1}(4)*params.w) + sin(deltaf)*(xi{k+1}(1) + params.lf*xi{k+1}(4))) <= e{k}(1),...
        params.rw*xi{k+1}(6) - (slipMax + 1)*(cos(deltaf)*(xi{k+1}(2) + xi{k+1}(4)*params.w) + sin(deltaf)*(xi{k+1}(1) + params.lf*xi{k+1}(4))) <= e{k}(2),...
        params.rw*xi{k+1}(7) - (xi{k+1}(2) - xi{k+1}(4)*params.w)*(slipMax + 1) <= e{k}(3),...
        params.rw*xi{k+1}(8) - (xi{k+1}(2) - xi{k+1}(4)*params.w)*(slipMax + 1) <= e{k}(4),...
        -(slipMax - 1)*(cos(deltaf)*(xi{k+1}(2) - xi{k+1}(4)*params.w) + sin(deltaf)*(xi{k+1}(1) + params.lf*xi{k+1}(4))) - params.rw*xi{k+1}(5) <= e{k}(5),...
        -(slipMax - 1)*(cos(deltaf)*(xi{k+1}(2) + xi{k+1}(4)*params.w) + sin(deltaf)*(xi{k+1}(1) + params.lf*xi{k+1}(4))) - params.rw*xi{k+1}(6) <= e{k}(6),...
        -(slipMax - 1)*(xi{k+1}(2) - xi{k+1}(4)*params.w) - params.rw*xi{k+1}(7) <= e{k}(7),...
        -(slipMax - 1)*(xi{k+1}(2) - xi{k+1}(4)*params.w) - params.rw*xi{k+1}(8) <= e{k}(8)];

    % tire slip angle constraints
    slipAngleConstraints = [        
        cos(deltaf)*(xi{k+1}(1) + params.lf*xi{k+1}(4)) - sin(deltaf)*(xi{k+1}(2) - xi{k+1}(4)*params.w) - slipAngleMax*(cos(deltaf)*(xi{k+1}(2) - xi{k+1}(4)*params.w) + sin(deltaf)*(xi{k+1}(1) + params.lf*xi{k+1}(4))) <= e{k}(9),...
        cos(deltaf)*(xi{k+1}(1) + params.lf*xi{k+1}(4)) - sin(deltaf)*(xi{k+1}(2) + xi{k+1}(4)*params.w) - slipAngleMax*(cos(deltaf)*(xi{k+1}(2) + xi{k+1}(4)*params.w) + sin(deltaf)*(xi{k+1}(1) + params.lf*xi{k+1}(4))) <= e{k}(10),...
        xi{k+1}(1) - params.lr*xi{k+1}(4) - slipAngleMax*(xi{k+1}(2) - xi{k+1}(4)*params.w) <= e{k}(11),...
        xi{k+1}(1) - params.lr*xi{k+1}(4) - slipAngleMax*(xi{k+1}(2) + xi{k+1}(4)*params.w) <= e{k}(12),...
        -slipAngleMax*(cos(deltaf)*(xi{k+1}(2) - xi{k+1}(4)*params.w) + sin(deltaf)*(xi{k+1}(1) + params.lf*xi{k+1}(4))) - cos(deltaf)*(xi{k+1}(1) + params.lf*xi{k+1}(4)) - sin(deltaf)*(xi{k+1}(2) - xi{k+1}(4)*params.w) <= e{k}(13),...
        -slipAngleMax*(cos(deltaf)*(xi{k+1}(2) + xi{k+1}(4)*params.w) + sin(deltaf)*(xi{k+1}(1) + params.lf*xi{k+1}(4))) - cos(deltaf)*(xi{k+1}(1) + params.lf*xi{k+1}(4)) - sin(deltaf)*(xi{k+1}(2) + xi{k+1}(4)*params.w) <= e{k}(14),...
        -slipAngleMax*(xi{k+1}(2) - xi{k+1}(4)*params.w) - xi{k+1}(1) - params.lr*xi{k+1}(4) <= e{k}(15),...
        -slipAngleMax*(xi{k+1}(2) + xi{k+1}(4)*params.w) - xi{k+1}(1) - params.lr*xi{k+1}(4) <= e{k}(16)];

    constraints = [constraints, slipConstraints, slipAngleConstraints];
    objective = objective + p * e{k}'*e{k};     
end
objective = objective + (C*xi{N+1}-xiref{N+1})'*QN*(C*xi{N+1}-xiref{N+1});

parameters_in = {xi{1},pastu,[xiref{:}],Ad,Bd,[xihat{:}],deltaf};
solutions_out = {[u{:}], [xi{:}]};            

controller = optimizer(constraints,objective,ops,parameters_in,solutions_out);

% Simulation length
nsim = 102;

% Initial and reference states
u0 = zeros(4,1);
vx0 = 10;
vx1 = 8;
xi0 = [0; vx0; 0; 0; 250/79*vx0*ones(4,1)];
xiref = repmat(C*xi0,1,nsim+N+1);  % constant forward driving
xiref = [zeros(1,nsim+N+1); 
         [vx0*ones(1,nsim/3),vx0+(vx1-vx0)/nsim*3*(1:nsim/3),vx1*ones(1,nsim/3+N+1)];
         zeros(2,nsim+N+1)];  % slow down from vx0 to vx1
deltaf = zeros(nsim,1);

% log data
Xis = nan(nx,nsim+1);
Us = nan(nu,nsim);
Xis(:,1) = xi0;

figure; hold on
Xi = xi0;
U = u0;
tic
for i = 1 : nsim
    % linearize around the current state and input
    [A,B] = twoTrackModelJacobian(Xi,U,deltaf(i));
    [Ad,Bd] = c2d(A,B,Ts);

    % calculate predictions with constant input
    [~,xihat] = ode45(@odeTwoTrackModel,0:Ts:N*Ts,[Xi;U;deltaf(i)]);
    xihat = xihat(:,1:8)';  % nx*(N+1)
    
    % solve the optimization problem
    inputs = {Xi,U,xiref(:,i:i+N),Ad,Bd,xihat,deltaf(i)};
    [solutions, diagnostics] = controller(inputs);
    U = solutions{1}(:,1);
    
    predictions = solutions{2};
    cla;
    stairs(i:i+N,xiref(2,i:i+N),'k')
    stairs(i:i+N,[Xi(2),predictions(2,2:end)],'b')
    stairs(1:i,Xis(2,1:i),'g');title('Prediction+reference+state')
    xlabel('k')
    ylabel('vx [m/s')
    legend('reference','prediction','trajectory')
    pause(0.01)
    
    % apply the inputs to the plant
    [~,Xi] = ode45(@odeTwoTrackModel,[0,Ts],[Xi;U;deltaf(i)]);
    Xi = Xi(end,1:8)';
    i
    
    % save data for plotting
    Xis(:,i+1) = Xi;
    Us(:,i) = U;
end
toc

% plot the results
t = 0:Ts:nsim*Ts;

% states and reference
figure  
plot(t,Xis')
hold on
plot(t,xiref(:,1:nsim+1)')
legend('v_y','v_x','\theta','der(\theta)',...
       '\omega_{fl}','\omega_{fr}','\omega_{rl}','\omega_{rr}',...
       'v_{y,ref}','v_{x,ref}','\theta_{ref}','der(\theta)_{ref}')
xlabel('t [s]')
ylabel('States')

% inputs
figure  
stairs(t,[u0 Us]')
legend('T_{fl}','T_{fr}','T_{rl}','T_{frr}')
xlabel('t [s]')
ylabel('Input torque [N]')

%% State constraints
        % state constraints
        % tire slip constraints
%         slipMax = 3.8/100;
%         slipConstraints = [    
%         params.rw*xi{k+1}(5) <= (slipMax + 1 + e{k}(1))*(cos(deltaf(i))*(xi{k+1}(2) - xi{k+1}(4)*params.w) + sin(deltaf(i))*(xi{k+1}(1) + params.lf*xi{k+1}(4))),...
%         params.rw*xi{k+1}(6) <= (slipMax + 1 + e{k}(2))*(cos(deltaf(i))*(xi{k+1}(2) + xi{k+1}(4)*params.w) + sin(deltaf(i))*(xi{k+1}(1) + params.lf*xi{k+1}(4))),...
%         params.rw*xi{k+1}(7) <= (xi{k+1}(2) - xi{k+1}(4)*params.w)*(slipMax + 1 + e{k}(3)),...
%         params.rw*xi{k+1}(8) <= (xi{k+1}(2) - xi{k+1}(4)*params.w)*(slipMax + 1 + e{k}(4)),...
%         -(e{k}(5) + slipMax - 1)*(cos(deltaf(i))*(xi{k+1}(2) - xi{k+1}(4)*params.w) + sin(deltaf(i))*(xi{k+1}(1) + params.lf*xi{k+1}(4))) <= params.rw*xi{k+1}(5),...
%         -(e{k}(6) + slipMax - 1)*(cos(deltaf(i))*(xi{k+1}(2) + xi{k+1}(4)*params.w) + sin(deltaf(i))*(xi{k+1}(1) + params.lf*xi{k+1}(4))) <= params.rw*xi{k+1}(6),...
%         -(xi{k+1}(2) - xi{k+1}(4)*params.w)*(e{k}(7) + slipMax - 1) <= params.rw*xi{k+1}(7),...
%         -(xi{k+1}(2) - xi{k+1}(4)*params.w)*(e{k}(8) + slipMax - 1) <= params.rw*xi{k+1}(8)];
%         
%         % tire slip angle constraints
%         slipAngleMax = deg2rad(3.4);
%         slipAngleConstraints = [        
%         cos(deltaf(i))*(xi{k+1}(1) + params.lf*xi{k+1}(4)) - sin(deltaf(i))*(xi{k+1}(2) - xi{k+1}(4)*params.w) <= (slipAngleMax + e{k}(9))*(cos(deltaf(i))*(xi{k+1}(2) - xi{k+1}(4)*params.w) + sin(deltaf(i))*(xi{k+1}(1) + params.lf*xi{k+1}(4))),...
%         cos(deltaf(i))*(xi{k+1}(1) + params.lf*xi{k+1}(4)) - sin(deltaf(i))*(xi{k+1}(2) + xi{k+1}(4)*params.w) <= (slipAngleMax + e{k}(10))*(cos(deltaf(i))*(xi{k+1}(2) + xi{k+1}(4)*params.w) + sin(deltaf(i))*(xi{k+1}(1) + params.lf*xi{k+1}(4))),...
%         xi{k+1}(1) - params.lr*xi{k+1}(4) <= (slipAngleMax + e{k}(11))*(xi{k+1}(2) - xi{k+1}(4)*params.w),...
%         xi{k+1}(1) - params.lr*xi{k+1}(4) <= (slipAngleMax + e{k}(12))*(xi{k+1}(2) + xi{k+1}(4)*params.w),...
%         -(slipAngleMax + e{k}(13))*(cos(deltaf(i))*(xi{k+1}(2) - xi{k+1}(4)*params.w) + sin(deltaf(i))*(xi{k+1}(1) + params.lf*xi{k+1}(4))) <= cos(deltaf(i))*(xi{k+1}(1) + params.lf*xi{k+1}(4)) - sin(deltaf(i))*(xi{k+1}(2) - xi{k+1}(4)*params.w),...
%         -(slipAngleMax + e{k}(14))*(cos(deltaf(i))*(xi{k+1}(2) + xi{k+1}(4)*params.w) + sin(deltaf(i))*(xi{k+1}(1) + params.lf*xi{k+1}(4))) <= cos(deltaf(i))*(xi{k+1}(1) + params.lf*xi{k+1}(4)) - sin(deltaf(i))*(xi{k+1}(2) + xi{k+1}(4)*params.w),...
%         -(slipAngleMax + e{k}(15))*(xi{k+1}(2) - xi{k+1}(4)*params.w) <= xi{k+1}(1) - params.lr*xi{k+1}(4),...
%         -(slipAngleMax + e{k}(16))*(xi{k+1}(2) + xi{k+1}(4)*params.w) <= xi{k+1}(1) - params.lr*xi{k+1}(4)];
%     
%         constraints = [constraints, slipConstraints, slipAngleConstraints];
%         objective = objective + p * e{k}'*e{k};
%% Pi-groups
%{
Mx = diag([sqrt(Cfx*l/m), sqrt(Cfx*l/m), 1, sqrt(Cfx/m/l)*ones(1,5)]);
Mxinv = Mx\eye(size(A,1));
Mu = l*Cfx*eye(4);
Muinv = Mu\eye(size(B,2));
that = sqrt(m*l/Cfx);  % time conversion factor t=that*tp

Ap = that*Mxinv*DfDx*Mx;
Bp = that*Mxinv*DfDu*Mu;
Cp = C*Mx;

Tsp = Ts/that;  % normalize time units
% [Apd, Bpd] = c2d(Ap,Bp,Tsp);

% from c2d.m
% [m,n] = size(Ap); %#ok<ASGLU>
% [m,nb] = size(Bp); %#ok<ASGLU>
% s = expm([[Ap Bp]*Tsp; zeros(nb,n+nb)]);
% Apd = s(1:n,1:n);
%}
% Bpd = s(1:n,n+1:n+nb);
