%% Initialize the simulation and create an MPC controller using YALMIP
clear;yalmip('clear');clc;close all

load ../../ParamsFull.mat
load ../yalmipsimulink/Params.mat
m   = VEHICLE.MASS;
L   = VEHICLE.WHEEL_BASE;
Cfx = VEHICLE.SLIP_STIFF;
w   = VEHICLE.TRACK_FRONT;
lf  = VEHICLE.LF;
lr  = VEHICLE.LR;
rw  = VEHICLE.WHEEL_RADIUS;
steering_ratio  = VEHICLE.STEERING_RATIO;
max_steering_rate = VEHICLE.MAX_STEERING_RATE;

% simulation parameters
Ts = 0.05;                  %[s] sampling time
tsim = 10;                  %[s] simulation time
nsim = round(tsim/Ts)+1;    %[-] simulation length
t = 0:Ts:tsim;              % time point vector

% initial state and input
SIM.VX_INIT = 20;  % m/s
SIM.U_INIT = 25;  % Nm
SIM.OMEGA_INIT = 250/79 * SIM.VX_INIT;  % rad/s

vx0 = SIM.VX_INIT;
xi0 = [vx0; 0; 0; 250/79*vx0*ones(4,1)];
u0 = SIM.U_INIT*ones(4,1);

% manoeuvre setup
tstart = 3;
manoeuvre = 1;
switch manoeuvre
    case 1  % step steer        
        steer_goal = 45;  % target steering angle in degrees
        steer_ramp = [0 : rad2deg(VEHICLE.MAX_STEERING_RATE)*Ts : steer_goal, steer_goal];  % limit the rate
        steer = [zeros(1,sum(t<tstart))  steer_ramp  steer_goal*ones(1,nsim-sum(t<tstart)-length(steer_ramp))];
        
    case 2  % sine steer
        steer_amplitude = 45;  % deg
        steer_frequency = 0.05;  % Hz
        steer = steer_amplitude*sin(2*pi*steer_frequency*(t - tstart));
        steer(t<tstart) = 0;
        
        steer_clip = 90;  % deg
        steer(steer>steer_clip) = steer_clip;
        steer(steer<-steer_clip) = -steer_clip;
        
    otherwise  % drive straight
        steer = zeros(size(t));
end
% change format for loading in Simulink
steer = [t' steer'];
vxref = [t' vx0*ones(size(t))'];
                  
% MPC data
nx = 7;  % no. of states
nu = 4;  % no. of inputs
C = [eye(3),zeros(3,4)];    % first three states are available
ny = size(C,1);             % no. of tracked states
N = 10;  % prediction horizon
M = 3;  % control horizon (blocking)

% Pi-groups
enablePi = 0;               % use Pi-groups?
Mxi = diag([sqrt(params.Cfx*params.l/params.m),...
           sqrt(params.Cfx*params.l/params.m),...
           sqrt(params.Cfx/params.m/params.l)*ones(1,5)]);
if ~enablePi; Mxi = eye(nx); end
Mxinv = Mxi\eye(nx);
Mu = params.l*params.Cfx*eye(nu);
if ~enablePi; Mu = eye(nu); end
Muinv = Mu\eye(nu);
that = sqrt(params.m*params.l/params.Cfx);  % time conversion factor t=that*tp
if ~enablePi; that = 1; end
Tsp = Ts/that;  % normalize time units

% create a data structure for the succesive linearization function
MPCData.nx = nx;
MPCData.nu = nu;
MPCData.ny = ny;
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

disp('Simulation init done')
%% create the controller

% select MPC setup
enable_constraints = [1;    % dynamics
                      0;    % state (slip, slip angle) -> soft constraints
                      1;    % input
                      1];   % input rate

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
Q = diag([1e2 0 1e4]);                  % tracked state weights
Q = C * Mxi * C' * Q * C * Mxi * C';    % transform to Pi-space
QN = Q;                                 % terminal cost weights
S = 1e-3*eye(nu);                       % input weights
S = Mu * S * Mu;                        % transform to Pi-space
R = 1e-2*eye(nu);                       % input rate weights
R = that^-2 * Mu * R * Mu;              % transform to Pi-space
p = 100;                                % slack variable weights

% slip and slip angle limits (for stability)
slipMax = VEHICLE.LIN_TIRE_KAPPA_MAX;  % conversion not needed (dimensionless)
slipMin = VEHICLE.LIN_TIRE_KAPPA_MIN;
alphaMax = VEHICLE.LIN_TIRE_ALPHA_MAX;  % conversion not needed (angle)
alphaMin = VEHICLE.LIN_TIRE_ALPHA_MIN;

% actuator limitations
torqueMax = VEHICLE.MAX_TORQUE;  % N
torqueMax = Muinv(1) * torqueMax;  % convert to Pi-space
torqueMin = Muinv(1) * -VEHICLE.MAX_TORQUE;

torqueRateMax = VEHICLE.MAX_TORQUE_RATE;  % Nm/s
torqueRateMax = Muinv(1) * that * torqueRateMax;  % convert to Pi-space
torqueRateMax = torqueRateMax * Tsp;  % max change in one time step
torqueRateMin = Muinv(1) * that * -VEHICLE.MAX_TORQUE_RATE * Tsp;

%% formulate the constraints
constraints = [];

% input constraint
if enable_constraints(3)
    constraints = [constraints, (torqueMin <= [u{:}] <= torqueMax):'Input limits'];
end

% input rate constraint
deltau = diff([pastu{1} u{:}],1,2);
if enable_constraints(4)
    constraints = [constraints, (torqueRateMin <= deltau <= torqueRateMax):'Input rate limits'];
end

% positive slack variables
if enable_constraints(2)
    constraints = [constraints, ([e{:}] >= 0):'Positive slack'];
end

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
        if k <= M && enable_constraints(4)
            objective = objective + deltau(:,k)' * R * deltau(:,k);
        end
    end
    
    % state constraints
    slipConstraints = [    
        (-cos(deltaf)*(slipMax + 1))*xi{k}(1) + (-sin(deltaf)*(slipMax + 1))*xi{k}(2) + (w*cos(deltaf) - lf*sin(deltaf))*(slipMax + 1)*xi{k}(3) + rw*xi{k}(4) - e{k}(1) <= 0,...
        (-cos(deltaf)*(slipMax + 1))*xi{k}(1) + (-sin(deltaf)*(slipMax + 1))*xi{k}(2) + (-(w*cos(deltaf) + lf*sin(deltaf))*(slipMax + 1))*xi{k}(3) + rw*xi{k}(5) - e{k}(2) <= 0,...
        (- slipMax - 1)*xi{k}(1) + w*(slipMax + 1)*xi{k}(3) + rw*xi{k}(6) - e{k}(3) <= 0,...
        (- slipMax - 1)*xi{k}(1) + (-w*(slipMax + 1))*xi{k}(3) + rw*xi{k}(7) - e{k}(4) <= 0,...
        cos(deltaf)*(slipMin + 1)*xi{k}(1) + sin(deltaf)*(slipMin + 1)*xi{k}(2) + (-(w*cos(deltaf) - lf*sin(deltaf))*(slipMin + 1))*xi{k}(3) + (-rw)*xi{k}(4) - e{k}(5) <= 0,...
        cos(deltaf)*(slipMin + 1)*xi{k}(1) + sin(deltaf)*(slipMin + 1)*xi{k}(2) + (w*cos(deltaf) + lf*sin(deltaf))*(slipMin + 1)*xi{k}(3) + (-rw)*xi{k}(5) - e{k}(6) <= 0,...
        (slipMin + 1)*xi{k}(1) + (-w*(slipMin + 1))*xi{k}(3) + (-rw)*xi{k}(6) - e{k}(7) <= 0,...
        (slipMin + 1)*xi{k}(1) + w*(slipMin + 1)*xi{k}(3) + (-rw)*xi{k}(7) - e{k}(8) <= 0];

    slipAngleConstraints = [
        (- sin(deltaf) - alphaMax*cos(deltaf))*xi{k}(1) + (cos(deltaf) - alphaMax*sin(deltaf))*xi{k}(2) + (alphaMax*(w*cos(deltaf) - lf*sin(deltaf)) + lf*cos(deltaf) + w*sin(deltaf))*xi{k}(3) - e{k}(9) <= 0,...
        (- sin(deltaf) - alphaMax*cos(deltaf))*xi{k}(1) + (cos(deltaf) - alphaMax*sin(deltaf))*xi{k}(2) + (lf*cos(deltaf) - alphaMax*(w*cos(deltaf) + lf*sin(deltaf)) - w*sin(deltaf))*xi{k}(3) - e{k}(10) <= 0,...
        (-alphaMax)*xi{k}(1) + xi{k}(2) + (alphaMax*w - lr)*xi{k}(3) - e{k}(11) <= 0,...
        (-alphaMax)*xi{k}(1) + xi{k}(2) + (- lr - alphaMax*w)*xi{k}(3) - e{k}(12) <= 0,...
        (sin(deltaf) + alphaMin*cos(deltaf))*xi{k}(1) + (alphaMin*sin(deltaf) - cos(deltaf))*xi{k}(2) + (- alphaMin*(w*cos(deltaf) - lf*sin(deltaf)) - lf*cos(deltaf) - w*sin(deltaf))*xi{k}(3) - e{k}(13) <= 0,...
        (sin(deltaf) + alphaMin*cos(deltaf))*xi{k}(1) + (alphaMin*sin(deltaf) - cos(deltaf))*xi{k}(2) + (alphaMin*(w*cos(deltaf) + lf*sin(deltaf)) - lf*cos(deltaf) + w*sin(deltaf))*xi{k}(3) - e{k}(14) <= 0,...
        alphaMin*xi{k}(1) - xi{k}(2) + (lr - alphaMin*w)*xi{k}(3) - e{k}(15) <= 0,...
        alphaMin*xi{k}(2) - xi{k}(2) + (lr + alphaMin*w)*xi{k}(3) - e{k}(16) <= 0];
    
    if enable_constraints(2)
        constraints = [constraints, slipConstraints:['Slip ' num2str(k)],...
                                    slipAngleConstraints:['Slip angle ' num2str(k)]];
    end
                            
    objective = objective + p * e{k}'*e{k};  % penalize slack variables     
end

parameters_in = {xi{1},[pastu{:}],[xiref{:}],[Ad{:}],[Bd{:}],[xihat{:}],deltaf};
solutions_out = {[u{:}], [xi{:}]};            

controller = optimizer(constraints,objective,ops,parameters_in,solutions_out);

save SavedController.mat controller

disp('Controller created')

%% test the controller
% load yalmipsimulink/Inputs.mat
% [solutions, errorcode, errortext] = controller(inputs);
% if errorcode ~= 0
%     error(errortext{1});
% else
%     disp('Controller test OK');
% end