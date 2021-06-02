%% Two-track model LTV MPC
clear;yalmip('clear');close all;clc
load params.mat  % vehicle parameters
nx = 8;  % no. of states
nu = 4;  % no. of inputs
C = [eye(4),zeros(4)];
ny = size(C,1);  % no. of tracked states
Ts = 0.05;  % sampling time
N = 10;  % prediction horizon
M = 3;  % control horizon (blocking)
nsim = 102;  % simulation length

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

% YALMIP data
Ad = sdpvar(nx*ones(N),nx*ones(N),'full');
Bd = sdpvar(nx*ones(N),nu*ones(N));
xihat = sdpvar(nx*ones(1,N+1),ones(1,N+1));  % predicted states from the previous step
u = sdpvar(nu*ones(1,M),ones(1,M));  % inputs
xi = sdpvar(nx*ones(1,N+1),ones(1,N+1));  % states
xiref = sdpvar(ny*ones(1,N+1),ones(1,N+1));  % reference
e = sdpvar(16*ones(1,N+1),ones(1,N+1));  % soft constraint variables
pastu = sdpvar(nu*ones(1,M),ones(1,M));  % predicted inputs from the previous step
deltaf = sdpvar(1,N+1);  % steering angle (for state constraints)
ops = sdpsettings('verbose',1,...
                  'solver','osqp',...
                  'osqp.eps_abs',1e-6,...
                  'osqp.eps_rel',1e-6,...
                  'debug',1,...
                  'cachesolvers',1);  % print output

% QP formulation
Q = diag([10 100e3 1e4 2e4]);  % tracked state weights
Q = C * Mxi * C' * Q * C * Mxi * C';  % transform to Pi-space
QN = Q;  % terminal cost weights
R = 1e-3*eye(nu);  % input weights
R = Mu * R * Mu;  % transform to Pi-space
p = 100;  % slack variable weights

slipMax = 3.8/100;  % conversion not needed
slipAngleMax = deg2rad(3.4);  % conversion not needed
torqueMax = 500;
torqueMax = Muinv(1) * torqueMax;  % convert to Pi-space
torqueRateMax = 50;
torqueRateMax = Muinv(1) * that * torqueRateMax;  % convert to Pi-space
constraints = [];
constraints = [constraints, -torqueMax <= [u{:}] <= torqueMax];  % input constraint
constraints = [constraints, -torqueRateMax <= diff([pastu{1} u{:}]) <= torqueRateMax];  % input rate constraint
constraints = [constraints, [e{:}] >= 0];  % slack variables are positive
objective = 0;
for k = 1:N+1
    if k == N+1
        objective = objective + (C*xi{k}-xiref{k})'*QN*(C*xi{k}-xiref{k});
    else
        if k > M
            objective = objective + u{M}'*R*u{M};
            constraints = [constraints, xi{k+1} == Ad{k}*xi{k} + Bd{k}*u{M} + xihat{k+1} - Ad{k}*xihat{k} - Bd{k}*pastu{M}];
        else
            objective = objective + u{k}'*R*u{k};
            constraints = [constraints, xi{k+1} == Ad{k}*xi{k} + Bd{k}*u{k} + xihat{k+1} - Ad{k}*xihat{k} - Bd{k}*pastu{k}];
        end
        objective = objective + (C*xi{k}-xiref{k})'*Q*(C*xi{k}-xiref{k});
    end
    
    % state constraints
    % tire slip constraints
    slipConstraints = [    
        params.rw*xi{k}(5) - (slipMax + 1)*(cos(deltaf(k))*(xi{k}(2) - xi{k}(4)*params.w) + sin(deltaf(k))*(xi{k}(1) + params.lf*xi{k}(4))) <= e{k}(1),...
        params.rw*xi{k}(6) - (slipMax + 1)*(cos(deltaf(k))*(xi{k}(2) + xi{k}(4)*params.w) + sin(deltaf(k))*(xi{k}(1) + params.lf*xi{k}(4))) <= e{k}(2),...
        params.rw*xi{k}(7) - (xi{k}(2) - xi{k}(4)*params.w)*(slipMax + 1) <= e{k}(3),...
        params.rw*xi{k}(8) - (xi{k}(2) - xi{k}(4)*params.w)*(slipMax + 1) <= e{k}(4),...
        -(slipMax - 1)*(cos(deltaf(k))*(xi{k}(2) - xi{k}(4)*params.w) + sin(deltaf(k))*(xi{k}(1) + params.lf*xi{k}(4))) - params.rw*xi{k}(5) <= e{k}(5),...
        -(slipMax - 1)*(cos(deltaf(k))*(xi{k}(2) + xi{k}(4)*params.w) + sin(deltaf(k))*(xi{k}(1) + params.lf*xi{k}(4))) - params.rw*xi{k}(6) <= e{k}(6),...
        -(slipMax - 1)*(xi{k}(2) - xi{k}(4)*params.w) - params.rw*xi{k}(7) <= e{k}(7),...
        -(slipMax - 1)*(xi{k}(2) - xi{k}(4)*params.w) - params.rw*xi{k}(8) <= e{k}(8)];

    % tire slip angle constraints
    slipAngleConstraints = [        
        cos(deltaf(k))*(xi{k}(1) + params.lf*xi{k}(4)) - sin(deltaf(k))*(xi{k}(2) - xi{k}(4)*params.w) - slipAngleMax*(cos(deltaf(k))*(xi{k}(2) - xi{k}(4)*params.w) + sin(deltaf(k))*(xi{k}(1) + params.lf*xi{k}(4))) <= e{k}(9),...
        cos(deltaf(k))*(xi{k}(1) + params.lf*xi{k}(4)) - sin(deltaf(k))*(xi{k}(2) + xi{k}(4)*params.w) - slipAngleMax*(cos(deltaf(k))*(xi{k}(2) + xi{k}(4)*params.w) + sin(deltaf(k))*(xi{k}(1) + params.lf*xi{k}(4))) <= e{k}(10),...
        xi{k}(1) - params.lr*xi{k}(4) - slipAngleMax*(xi{k}(2) - xi{k}(4)*params.w) <= e{k}(11),...
        xi{k}(1) - params.lr*xi{k}(4) - slipAngleMax*(xi{k}(2) + xi{k}(4)*params.w) <= e{k}(12),...
        -slipAngleMax*(cos(deltaf(k))*(xi{k}(2) - xi{k}(4)*params.w) + sin(deltaf(k))*(xi{k}(1) + params.lf*xi{k}(4))) - cos(deltaf(k))*(xi{k}(1) + params.lf*xi{k}(4)) - sin(deltaf(k))*(xi{k}(2) - xi{k}(4)*params.w) <= e{k}(13),...
        -slipAngleMax*(cos(deltaf(k))*(xi{k}(2) + xi{k}(4)*params.w) + sin(deltaf(k))*(xi{k}(1) + params.lf*xi{k}(4))) - cos(deltaf(k))*(xi{k}(1) + params.lf*xi{k}(4)) - sin(deltaf(k))*(xi{k}(2) + xi{k}(4)*params.w) <= e{k}(14),...
        -slipAngleMax*(xi{k}(2) - xi{k}(4)*params.w) - xi{k}(1) - params.lr*xi{k}(4) <= e{k}(15),...
        -slipAngleMax*(xi{k}(2) + xi{k}(4)*params.w) - xi{k}(1) - params.lr*xi{k}(4) <= e{k}(16)];

    constraints = [constraints, slipConstraints, slipAngleConstraints];
    objective = objective + p * e{k}'*e{k};     
end

parameters_in = {xi{1},[pastu{:}],[xiref{:}],[Ad{:}],[Bd{:}],[xihat{:}],deltaf};
solutions_out = {[u{:}], [xi{:}]};            

controller = optimizer(constraints,objective,ops,parameters_in,solutions_out);

% Initial and reference states
u0 = 3.2590*ones(4,1);
vx0 = 10;
vx1 = 8;
xi0 = [0; vx0; 0; 0; 250/79*vx0*ones(4,1)];
xiref = repmat(C*xi0,1,nsim+N+1);  % constant forward driving
xiref = [zeros(1,nsim+N+1); 
         [vx0*ones(1,nsim/3),vx0+(vx1-vx0)/nsim*3*(1:nsim/3),vx1*ones(1,nsim/3+N+1)];
         zeros(2,nsim+N+1)];  % slow down from vx0 to vx1
xiref = C * Mxinv * C' * xiref;  % convert to Pi-space
deltaf = zeros(1,nsim+N);

% linearization data from the previous step
oldXis = repmat(xi0,1,N+1);
oldUs = repmat(u0,1,M);
xihats = nan(nx,N+1);

% Jacobian matrices
As = cell(1,N);
%As = nan(nx,nx,N+1);
Bs = cell(1,N);
%Bs = nan(nx,nu,N+1);

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
    for j = 1:N
        if j > M
            [~,xihat] = ode45(@odeTwoTrackModel,[0,Ts],[xihat;oldUs(:,M);deltaf(i+j-1)]);
            [A,B] = twoTrackModelJacobian(oldXis(:,j),oldUs(:,M),deltaf(i+j-1));
        else
            [~,xihat] = ode45(@odeTwoTrackModel,[0,Ts],[xihat;oldUs(:,j);deltaf(i+j-1)]);
            [A,B] = twoTrackModelJacobian(oldXis(:,j),oldUs(:,j),deltaf(i+j-1));
        end
        xihat = xihat(end,1:8)';  % nx*(N+1)
        % transform to Pi-space
        xihats(:,j+1) = Mxinv * xihat;
        
        % formulate in Pi-space and discretize
        Ap = that * Mxinv * A * Mxi;
        Bp = that * Mxinv * B * Mu;
        [Ad,Bd] = c2d(Ap,Bp,Tsp);
        As{j} = Ad;
        Bs{j} = Bd;
        %As(:,:,j) = Ad;
        %Bs(:,:,j) = Bd;        
    end
    oldUs = Muinv * oldUs;
    
    Xi = Mxinv * Xi;
    U = Muinv * U;
    
    % solve the optimization problem
    inputs = {Xi,oldUs,xiref(:,i:i+N),[As{:}],[Bs{:}],xihats,deltaf(i:i+N)};
    [solutions, diagnostics] = controller(inputs);
    U = solutions{1}(:,1);  % take the first optimal input
    
    if diagnostics == 1
        error('The problem is infeasible');
    end
    
    % "live" plot
    predictions = solutions{2};
    cla;
    stairs(i:i+N,Mxi(2,2)*xiref(2,i:i+N),'k')
    stairs(i:i+N,Mxi(2,2)*[Xi(2),predictions(2,2:end)],'b')
    stairs(1:i,Xis(2,1:i),'g');title('Prediction+reference+state')
    xlabel('k')
    ylabel('vx [m/s')
    legend('reference','prediction','trajectory')
    pause(0.01)
    
    % apply the inputs to the plant
    U = Mu * U
    Xi = Mxi * Xi;
    [~,Xi] = ode45(@odeSimulationModel,[0,Ts],[Xi;U;deltaf(i);ay;ax]);
    Xi = Xi(end,1:8)'
    i
    
    % calculate the acceleration
    f = odeSimulationModel(0,[Xi;U;deltaf(i);ay;ax]);
    ay = f(1); ax = f(2);
    
    % save optimization results for the next iteration
    oldXis = [Xi, Mxi * predictions(:,2:end)];
    oldUs = Mu * solutions{1};
    
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
