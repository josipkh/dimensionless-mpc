%% Initialize the simulation and create an MPC controller using YALMIP
clear;yalmip('clear');clc;close all

load ../ParamsFull.mat
load yalmipsimulink/Params.mat

% simulation parameters
Ts = 0.05;                  %[s] sampling time
tsim = 100;                 %[s] simulation time
nsim = round(tsim/Ts)+1;    %[-] simulation length
t = 0:Ts:tsim;              % time point vector

% MPC data
nx = 8;  % no. of states
nu = 4;  % no. of inputs
C = [eye(4),zeros(4)];
ny = size(C,1);  % no. of tracked states
N = 10;  % prediction horizon
M = 3;  % control horizon (blocking)

% initial state and input
vx0 = 10;  %[m/s]
xi0 = [0; vx0; 0; 0; 250/79*vx0*ones(4,1)];
u0 = 3.2590*ones(4,1);  %[Nm]

% reference generation
reference = 3;
switch reference
    case 1
        % constant forward driving
        u0 = zeros(4,1);
        deltaf = zeros(1,nsim);
        xiref = repmat(C*xi0,1,nsim+N+1);
    case 2
        % slow down from vx0 to vx1
        u0 = zeros(4,1);
        deltaf = zeros(1,nsim);
        vx1 = 8;
        xiref = [zeros(1,nsim+N+1); 
                 [vx0*ones(1,floor(nsim/3)),...
                  vx0+(vx1-vx0)/round(nsim/3)*(0:floor(nsim/3)),...
                  vx1*ones(1,nsim-2*floor(nsim/3)+N)];
                 zeros(2,nsim+N+1)]; 
    case 3
        % step steer
        steerStep = 5;  % steering angle in degrees
        deltaf = [zeros(1,floor(nsim/5)) deg2rad(steerStep)*ones(1,ceil(nsim*4/5))];  % step steer
        thetadref = vx0/params.l * tan(deltaf);  % eq. 23 from (FER,2019)
        thetaref = Ts * cumsum(thetadref);  % integration to get theta
        xiref = [zeros(1,nsim+N+1);
                 vx0*ones(1,nsim+N+1);
                 thetaref, thetaref(end)*ones(1,N+1);
                 thetadref, thetadref(end)*ones(1,N+1)];
    otherwise
        error('undefined reference')
end

% create time series for the reference
deltafs = [t',deltaf'];
xirefs = [t',nan(length(t),ny*(N+1))];
for i = 1:length(t)
    xiref_temp = xiref(:,i:i+N);
    xirefs(i,2:end) = xiref_temp(:);
end
% xirefs = [t',xiref(:,1:nsim)'];
disp('Reference created')
%==================================================================

% Pi-groups
enablePi = 1;
Mxi = diag([sqrt(params.Cfx*params.l/params.m),...
           sqrt(params.Cfx*params.l/params.m),...
           1,...
           sqrt(params.Cfx/params.m/params.l)*ones(1,5)]);
if ~enablePi; Mxi = eye(8); end
Mxinv = Mxi\eye(nx);
Mu = params.l*params.Cfx*eye(nu);
if ~enablePi; Mu = eye(4); end
Muinv = Mu\eye(nu);
that = sqrt(params.m*params.l/params.Cfx);  % time conversion factor t=that*tp
if ~enablePi; that = 1; end
Tsp = Ts/that;  % normalize time units

% create a data structure for the succesive linearization function
MPCData.nx = nx;
MPCData.nu = nu;
MPCData.N = N;
MPCData.M = M;
MPCData.Mxi = Mxi;
MPCData.Mxinv = Mxinv;
MPCData.Mu = Mu;
MPCData.Muinv = Muinv;
MPCData.that = that;
MPCData.Ts = Ts;
MPCData.Tsp = Tsp;
MPCData.params = params;

% convert reference to Pi-space
% xirefs = C * Mxinv * C' * xiref;
% xirefs = [t',xirefs(:,1:nsim)'];

%% create the controller

% YALMIP data
Ad = sdpvar(nx*ones(N),nx*ones(N),'full');
Bd = sdpvar(nx*ones(N),nu*ones(N));
xihat = sdpvar(nx*ones(1,N+1),ones(1,N+1));  % predicted states from the previous step
u = sdpvar(nu*ones(1,M),ones(1,M));  % inputs
xi = sdpvar(nx*ones(1,N+1),ones(1,N+1));  % states
xiref = sdpvar(ny*ones(1,N+1),ones(1,N+1));  % reference
e = sdpvar(16*ones(1,N+1),ones(1,N+1));  % soft constraint variables
pastu = sdpvar(nu*ones(1,M),ones(1,M));  % predicted inputs from the previous step
%deltaf = sdpvar(1,N+1);  % steering angle (for state constraints)
deltaf = sdpvar(1,1);  % assume that steering angle is constant during the horizon
ops = sdpsettings('verbose',0,...
                  'solver','osqp',...
                  'osqp.eps_abs',1e-6,...
                  'osqp.eps_rel',1e-6,...
                  'debug',0,...
                  'cachesolvers',1);  % print output
% ops.savesolveroutput = 1;
% ops.savesolverinput = 1;
% ops.usex0 = 1;
% ops = sdpsettings('solver','quadprog');

% QP formulation
Q = diag([10 100 1e4 2e4]);  % tracked state weights
Q = C * Mxi * C' * Q * C * Mxi * C';  % transform to Pi-space
QN = Q;  % terminal cost weights
S = 1e-3*eye(nu);  % input weights
S = Mu * S * Mu;  % transform to Pi-space
R = 1e-2*eye(nu);  % input rate weights
R = that^-2 * Mu * R * Mu;  % transform to Pi-space
p = 100;  % slack variable weights

% slip and slip angle limits (for stability)
slipMax = params.slipMax;  % conversion not needed (dimensionless)
slipAngleMax = params.alphaMax;  % conversion not needed (angle)
% actuator limitations
torqueMax = params.Tmax;  % N
torqueMax = Muinv(1) * torqueMax;  % convert to Pi-space
torqueRateMax = params.dTmax;  % Nm
torqueRateMax = Muinv(1) * that * torqueRateMax;  % convert to Pi-space
constraints = [];
constraints = [constraints, (-torqueMax <= [u{:}] <= torqueMax):'Input limits'];
deltau = diff([pastu{1} u{:}],1,2);
constraints = [constraints, (-torqueRateMax <= deltau <= torqueRateMax):'Input rate limits'];
constraints = [constraints, ([e{:}] >= 0):'Positive slack'];
objective = 0;
for k = 1:N+1
    % system dynamics and the objective function
    if k == N+1
        objective = objective + (C*xi{k}-xiref{k})'*QN*(C*xi{k}-xiref{k});
    else
        if k > M
            objective = objective + u{M}'*S*u{M};
            constraints = [constraints, (xi{k+1} == Ad{k}*xi{k} + Bd{k}*u{M} + xihat{k+1} - Ad{k}*xihat{k} - Bd{k}*pastu{M}):['Dynamics ' num2str(k)]];
        else
            objective = objective + u{k}'*S*u{k};
            constraints = [constraints, (xi{k+1} == Ad{k}*xi{k} + Bd{k}*u{k} + xihat{k+1} - Ad{k}*xihat{k} - Bd{k}*pastu{k}):['Dynamics ' num2str(k)]];
        end
        objective = objective + (C*xi{k}-xiref{k})'*Q*(C*xi{k}-xiref{k});
        if k < M
            objective = objective + deltau(:,k)' * R * deltau(:,k);
        end
    end
    
    % state constraints
    % tire slip constraints
    slipConstraints = [    
        params.rw*xi{k}(5) - (slipMax + 1)*(cos(deltaf)*(xi{k}(2) - xi{k}(4)*params.w) + sin(deltaf)*(xi{k}(1) + params.lf*xi{k}(4))) <= e{k}(1),...
        params.rw*xi{k}(6) - (slipMax + 1)*(cos(deltaf)*(xi{k}(2) + xi{k}(4)*params.w) + sin(deltaf)*(xi{k}(1) + params.lf*xi{k}(4))) <= e{k}(2),...
        params.rw*xi{k}(7) - (xi{k}(2) - xi{k}(4)*params.w)*(slipMax + 1) <= e{k}(3),...
        params.rw*xi{k}(8) - (xi{k}(2) - xi{k}(4)*params.w)*(slipMax + 1) <= e{k}(4),...
        -(slipMax - 1)*(cos(deltaf)*(xi{k}(2) - xi{k}(4)*params.w) + sin(deltaf)*(xi{k}(1) + params.lf*xi{k}(4))) - params.rw*xi{k}(5) <= e{k}(5),...
        -(slipMax - 1)*(cos(deltaf)*(xi{k}(2) + xi{k}(4)*params.w) + sin(deltaf)*(xi{k}(1) + params.lf*xi{k}(4))) - params.rw*xi{k}(6) <= e{k}(6),...
        -(slipMax - 1)*(xi{k}(2) - xi{k}(4)*params.w) - params.rw*xi{k}(7) <= e{k}(7),...
        -(slipMax - 1)*(xi{k}(2) - xi{k}(4)*params.w) - params.rw*xi{k}(8) <= e{k}(8)];

    % tire slip angle constraints
    slipAngleConstraints = [        
        cos(deltaf)*(xi{k}(1) + params.lf*xi{k}(4)) - sin(deltaf)*(xi{k}(2) - xi{k}(4)*params.w) - slipAngleMax*(cos(deltaf)*(xi{k}(2) - xi{k}(4)*params.w) + sin(deltaf)*(xi{k}(1) + params.lf*xi{k}(4))) <= e{k}(9),...
        cos(deltaf)*(xi{k}(1) + params.lf*xi{k}(4)) - sin(deltaf)*(xi{k}(2) + xi{k}(4)*params.w) - slipAngleMax*(cos(deltaf)*(xi{k}(2) + xi{k}(4)*params.w) + sin(deltaf)*(xi{k}(1) + params.lf*xi{k}(4))) <= e{k}(10),...
        xi{k}(1) - params.lr*xi{k}(4) - slipAngleMax*(xi{k}(2) - xi{k}(4)*params.w) <= e{k}(11),...
        xi{k}(1) - params.lr*xi{k}(4) - slipAngleMax*(xi{k}(2) + xi{k}(4)*params.w) <= e{k}(12),...
        -slipAngleMax*(cos(deltaf)*(xi{k}(2) - xi{k}(4)*params.w) + sin(deltaf)*(xi{k}(1) + params.lf*xi{k}(4))) - cos(deltaf)*(xi{k}(1) + params.lf*xi{k}(4)) - sin(deltaf)*(xi{k}(2) - xi{k}(4)*params.w) <= e{k}(13),...
        -slipAngleMax*(cos(deltaf)*(xi{k}(2) + xi{k}(4)*params.w) + sin(deltaf)*(xi{k}(1) + params.lf*xi{k}(4))) - cos(deltaf)*(xi{k}(1) + params.lf*xi{k}(4)) - sin(deltaf)*(xi{k}(2) + xi{k}(4)*params.w) <= e{k}(14),...
        -slipAngleMax*(xi{k}(2) - xi{k}(4)*params.w) - xi{k}(1) - params.lr*xi{k}(4) <= e{k}(15),...
        -slipAngleMax*(xi{k}(2) + xi{k}(4)*params.w) - xi{k}(1) - params.lr*xi{k}(4) <= e{k}(16)];

    constraints = [constraints, slipConstraints:['Slip ' num2str(k)],...
                                slipAngleConstraints:['Slip angle ' num2str(k)]];
                            
    objective = objective + p * e{k}'*e{k};  % penalize slack variables     
end

parameters_in = {xi{1},[pastu{:}],[xiref{:}],[Ad{:}],[Bd{:}],[xihat{:}],deltaf};
solutions_out = {[u{:}], [xi{:}]};            

controller = optimizer(constraints,objective,ops,parameters_in,solutions_out);

save SavedController.mat controller

disp('Controller created')
%% test the controller
load yalmipsimulink/Inputs.mat
[solutions, errorcode, errortext] = controller(inputs);
if errorcode ~= 0
    error(errortext{1});
else
    disp('Controller test OK');
end