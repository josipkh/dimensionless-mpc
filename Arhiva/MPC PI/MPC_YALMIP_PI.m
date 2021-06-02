%% Two-track model LTV MPC
% TODO: control horizon, soft constraints, state constraints
clear;yalmip('clear');close all;clc
load params.mat  % vehicle parameters
nx = 8;  % no. of states
nu = 4;  % no. of inputs
Ts = 0.05;  % sampling time
N = 10;  % horizon length
M = 3;  % control horizon (blocking)

% Pi-groups
Mxi = diag([sqrt(params.Cfx*params.l/params.m),...
           sqrt(params.Cfx*params.l/params.m),...
           1,...
           sqrt(params.Cfx/params.m/params.l)*ones(1,5)]);
Mxinv = Mxi\eye(nx);
Mu = params.l*params.Cfx*eye(nu);
Muinv = Mu\eye(nu);
that = sqrt(params.m*params.l/params.Cfx);  % time conversion factor t=that*tp
Tsp = Ts/that;  % normalize time units

C = [eye(4),zeros(4)];
Cp = C * Mxi;
ny = size(C,1);  % no. of tracked states

Q = diag([10 100e3 1e4 2e4]);  % tracked state weights
Q = C * Mxi * C' * Q * C * Mxi * C';
QN = Q;  % terminal cost weights
R = 1e-3*eye(nu);  % input weights
R = Mu * R * Mu;

% YALMIP data
u = sdpvar(nu*ones(1,N),ones(1,N));  % inputs
xi = sdpvar(nx*ones(1,N+1),ones(1,N+1));  % states
%e = sdpvar(nx*ones(1,N+1),ones(1,N+1));  % soft constraint variables
pastu = sdpvar(4,1);  % for input rate limits
ops = sdpsettings('verbose',1,...
                  'solver','osqp',...
                  'osqp.eps_abs',1e-6,...
                  'osqp.eps_rel',1e-6,...
                  'debug',1,...
                  'cachesolvers',1);  % print output

% Simulation length
nsim = 30;

% Initial and reference states
vx0 = 10;  % m/s
vx1 = 8;
xi0 = [0; vx0; 0; 0; 250/79*vx0*ones(4,1)];
xiref = repmat(C*xi0,1,nsim+N);  % constant forward driving
% xiref = [zeros(1,nsim+N); 
%          [vx0*ones(1,nsim/3),vx0+(vx1-vx0)/nsim*3*(1:nsim/3),vx1*ones(1,nsim/3+N)];
%          zeros(2,nsim+N)];
% xiref = [zeros(1,nsim+N); 
%          [vx0*ones(1,20),vx1*ones(1,nsim+N-20)];
%          zeros(2,nsim+N)];
u0 = 3.2590*ones(4,1);
deltaf = zeros(nsim+N,1);

xiref = C * Mxinv * C' * xiref;

% linearization data from the previous step
oldXis = repmat(xi0,1,N+1);
oldUs = repmat(u0,1,N+1);
xihats = nan(nx,N+1);

% Jacobian matrices
As = cell(1,N+1);
Bs = cell(1,N+1);

% log data
Xis = nan(nx,nsim+1);
Us = nan(nu,nsim);
Xis(:,1) = xi0;

figure; hold on
Xi = xi0;
U = u0;
ay = 0;
ax = 0;

tic
for i = 1 : nsim    
    % linearize around the predicted states and inputs
    xihat = Xi;
    xihats(:,1) = Mxinv * Xi;
    for j = 1:N+1
        [~,xihat] = ode45(@odeTwoTrackModel,[0,Ts],[xihat;oldUs(:,j);deltaf(i+j-1)]);
        xihat = xihat(end,1:8)';
        % transform to Pi-space
        xihats(:,j+1) = Mxinv * xihat;
        
        [A,B] = twoTrackModelJacobian(oldXis(:,j),oldUs(:,j),deltaf(i+j-1));
        % formulate in Pi-space and discretize
        Ap = that * Mxinv * A * Mxi;
        Bp = that * Mxinv * B * Mu;
        [Ad,Bd] = c2d(Ap,Bp,Tsp);
        As{j} = Ad;
        Bs{j} = Bd;
    end
    oldUs = Muinv * oldUs;
    
    % constraints and cost
    Xi = Mxinv * Xi;
    U = Muinv * U;
    constraints = [xi{1} == Xi];  % initial state constraint
    %constraints = [constraints, pastu == U];
    %constraints = [constraints, -500 <= [u{:}] <= 500];  % input constraint
    %constraints = [constraints, -50 <= diff([pastu u{:}]) <= 50];  % input rate constraint    
    objective = 0;
    for k = 1:N
        % tracking cost
        objective = objective + (C*xi{k}-xiref(:,i+k-1))'*Q*(C*xi{k}-xiref(:,i+k-1));
        objective = objective + u{k}'*R*u{k};
        
        % calculate linearization offset
        d = xihats(:,k+1) - As{k}*xihats(:,k) - Bs{k}*oldUs(:,k);
        
        % system dynamics
        constraints = [constraints, xi{k+1} == As{k}*xi{k} + Bs{k}*u{k} + d];
    end
    % terminal cost
    objective = objective + (C*xi{N+1}-xiref(:,i+N))'*QN*(C*xi{N+1}-xiref(:,i+N));
    
    % solve the optimization problem
    diagnostics = optimize(constraints,objective,ops);
    U = value(u{1});
    
    predictions = value([xi{:}]);
    cla;
    stairs(i:i+N,Mxi(2,2)*xiref(2,i:i+N),'k')
    stairs(i:i+N,Mxi(2,2)*predictions(2,:),'b')
    stairs(1:i,Xis(2,1:i),'g');title('Prediction+reference+state')
    xlabel('k')
    ylabel('vx [m/s')
    legend('reference','prediction','trajectory')
    pause(0.01)
    
    % Apply the nu inputs to the plant
    U = Mu * U
    Xi = Mxi * Xi;
    
    [~,Xi] = ode45(@odeSimulationModel,[0,Ts],[Xi;U;deltaf(i);ay;ax]);
    Xi = Xi(end,1:8)'
    i
    f = odeSimulationModel(0,[Xi;U;deltaf(i);ay;ax]);
    ay = f(1); ax = f(2);
    oldXis = [Xi, Mxi * predictions(:,2:end)];
    oldUs = Mu * [value([u{:}]), value(u{N})];
    
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
plot(t,(Mxi(1:4,1:4)*xiref(:,1:nsim+1))')
legend('v_y','v_x','\theta','der(\theta)',...
       '\omega_{fl}','\omega_{fr}','\omega_{rl}','\omega_{rr}',...
       'v_{y,ref}','v_{x,ref}','\theta_{ref}','der(\theta)_{ref}')
xlabel('t [s]')
ylabel('States')

% inputs
figure  
plot(t,[u0 Us])
legend('T_{fl}','T_{fr}','T_{rl}','T_{frr}')
xlabel('t [s]')
ylabel('Input torque [N]')
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
%%
        % state constraints
        % tire slip constraints
%         slipMax = 3.8/100;
%         slipConstraints = [    
%         params.rw*xi{k+1}(5) <= (slipMax + 1)*(cos(deltaf(i))*(xi{k+1}(2) - xi{k+1}(4)*params.w) + sin(deltaf(i))*(xi{k+1}(1) + params.lf*xi{k+1}(4))),...
%         params.rw*xi{k+1}(6) <= (slipMax + 1)*(cos(deltaf(i))*(xi{k+1}(2) + xi{k+1}(4)*params.w) + sin(deltaf(i))*(xi{k+1}(1) + params.lf*xi{k+1}(4))),...
%         params.rw*xi{k+1}(7) <= (xi{k+1}(2) - xi{k+1}(4)*params.w)*(slipMax + 1),...
%         params.rw*xi{k+1}(8) <= (xi{k+1}(2) - xi{k+1}(4)*params.w)*(slipMax + 1),...
%         -(slipMax - 1)*(cos(deltaf(i))*(xi{k+1}(2) - xi{k+1}(4)*params.w) + sin(deltaf(i))*(xi{k+1}(1) + params.lf*xi{k+1}(4))) <= params.rw*xi{k+1}(5),...
%         -(slipMax - 1)*(cos(deltaf(i))*(xi{k+1}(2) + xi{k+1}(4)*params.w) + sin(deltaf(i))*(xi{k+1}(1) + params.lf*xi{k+1}(4))) <= params.rw*xi{k+1}(6),...
%         -(xi{k+1}(2) - xi{k+1}(4)*params.w)*(slipMax - 1) <= params.rw*xi{k+1}(7),...
%         -(xi{k+1}(2) - xi{k+1}(4)*params.w)*(slipMax - 1) <= params.rw*xi{k+1}(8)];
%         
%         % tire slip angle constraints
%         slipAngleMax = deg2rad(3.4);
%         slipAngleConstraints = [        
%         cos(deltaf(i))*(xi{k+1}(1) + params.lf*xi{k+1}(4)) - sin(deltaf(i))*(xi{k+1}(2) - xi{k+1}(4)*params.w) <= slipAngleMax*(cos(deltaf(i))*(xi{k+1}(2) - xi{k+1}(4)*params.w) + sin(deltaf(i))*(xi{k+1}(1) + params.lf*xi{k+1}(4))),...
%         cos(deltaf(i))*(xi{k+1}(1) + params.lf*xi{k+1}(4)) - sin(deltaf(i))*(xi{k+1}(2) + xi{k+1}(4)*params.w) <= slipAngleMax*(cos(deltaf(i))*(xi{k+1}(2) + xi{k+1}(4)*params.w) + sin(deltaf(i))*(xi{k+1}(1) + params.lf*xi{k+1}(4))),...
%         xi{k+1}(1) - params.lr*xi{k+1}(4) <= slipAngleMax*(xi{k+1}(2) - xi{k+1}(4)*params.w),...
%         xi{k+1}(1) - params.lr*xi{k+1}(4) <= slipAngleMax*(xi{k+1}(2) + xi{k+1}(4)*params.w),...
%         -slipAngleMax*(cos(deltaf(i))*(xi{k+1}(2) - xi{k+1}(4)*params.w) + sin(deltaf(i))*(xi{k+1}(1) + params.lf*xi{k+1}(4))) <= cos(deltaf(i))*(xi{k+1}(1) + params.lf*xi{k+1}(4)) - sin(deltaf(i))*(xi{k+1}(2) - xi{k+1}(4)*params.w),...
%         -slipAngleMax*(cos(deltaf(i))*(xi{k+1}(2) + xi{k+1}(4)*params.w) + sin(deltaf(i))*(xi{k+1}(1) + params.lf*xi{k+1}(4))) <= cos(deltaf(i))*(xi{k+1}(1) + params.lf*xi{k+1}(4)) - sin(deltaf(i))*(xi{k+1}(2) + xi{k+1}(4)*params.w),...
%         -slipAngleMax*(xi{k+1}(2) - xi{k+1}(4)*params.w) <= xi{k+1}(1) - params.lr*xi{k+1}(4),...
%         -slipAngleMax*(xi{k+1}(2) + xi{k+1}(4)*params.w) <= xi{k+1}(1) - params.lr*xi{k+1}(4)];
%     
%         constraints = [constraints, slipConstraints, slipAngleConstraints];        