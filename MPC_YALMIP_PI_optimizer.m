%% Two-track model LTV MPC
%% Problem setup
clear;yalmip('clear');close all;clc

% vehicle parameters
load ParamsFull.mat
m   = VEHICLE.MASS;
L   = VEHICLE.WHEEL_BASE;
Cfx = VEHICLE.SLIP_STIFF;
w   = VEHICLE.TRACK_FRONT;
lf  = VEHICLE.LF;
lr  = VEHICLE.LR;
rw  = VEHICLE.WHEEL_RADIUS;
steering_ratio  = VEHICLE.STEERING_RATIO;
max_steering_rate = VEHICLE.MAX_STEERING_RATE;

% system dimensions
nx = 7;                     % no. of states (vx, vy, dtheta, wheel speeds)
nu = 4;                     % no. of inputs (wheel torques)
neps = 16;                  % no. of slack variables
C = [eye(3),zeros(3,4)];    % first three states are available
ny = size(C,1);             % no. of tracked states

% simulation setup
Ts = 0.05;                  % sampling time [s]
Tsim = 10;                   % simulation time [s]
nsim = Tsim / Ts;           % no. of simulation steps

% MPC parameters
N = 10;                     % prediction horizon
M = 3;                      % control horizon (blocking)

enablePi = 0;               % use Pi-groups?
enable_constraints = [1;    % dynamics
                      0;    % state (slip, slip angle) -> soft constraints
                      1;    % input
                      0];   % input rate

% Pi-groups
Mxi = diag([sqrt(Cfx*L/m),...
           sqrt(Cfx*L/m),...
           sqrt(Cfx/m/L)*ones(1,5)]);
if ~enablePi; Mxi = eye(nx); end
Mxinv = Mxi\eye(nx);
Mu = L*Cfx * eye(nu);
if ~enablePi; Mu = eye(nu); end
Muinv = Mu\eye(nu);
that = sqrt(m*L/Cfx);  % time conversion factor t=that*tp
if ~enablePi; that = 1; end
Tsp = Ts/that;  % normalize time units

% objective function matrices
Q = diag([1e6 0 1e8]);                  % tracked state weights
Q = C * Mxi * C' * Q * C * Mxi * C';    % transform to Pi-space
QN = Q;                                 % terminal cost weights
S = 1e-3*eye(nu);                       % input weights
S = Mu * S * Mu;                        % transform to Pi-space
R = 1e-2*eye(nu);                       % input rate weights
R = that^-2 * Mu * R * Mu;              % transform to Pi-space
p = 100;                                % slack variable weights

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
ops = sdpsettings('verbose',0,... % print output? (0-2)
                  'solver','osqp',...
                  'osqp.eps_abs',1e-6,...
                  'osqp.eps_rel',1e-6,...
                  'debug',1,...
                  'cachesolvers',1);  
ops.savesolveroutput = 1;
ops.savesolverinput = 1;

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

%% initial state and reference

% initial state and input
vx0 = 10;
xi0 = [vx0; 0; 0; 250/79*vx0*ones(4,1)];
u0 = 3.2590*ones(4,1);

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
        xiref = [[vx0*ones(1,floor(nsim/3)),...
                  vx0+(vx1-vx0)/round(nsim/3)*(0:floor(nsim/3)),...
                  vx1*ones(1,nsim-2*floor(nsim/3)+N)];
                 zeros(2,nsim+N+1)]; 
    case 3
        % step steer (rate limited)
        steerGoal = 45;  % target steering angle in degrees
        steerGoal = deg2rad(steerGoal);  % convert to radians
        steerRamp = [0 : max_steering_rate*Ts : steerGoal, steerGoal];  % limit the rate
        deltaf = [zeros(1,floor(nsim/5)) steerRamp steerGoal*ones(1,nsim-floor(nsim/5)-length(steerRamp))];
        thetadref = vx0/L * tan(deltaf / steering_ratio);  % eq. 23 from (FER,2019)
        xiref = [vx0*ones(1,nsim+N+1);
                 zeros(1,nsim+N+1);
                 thetadref, thetadref(end)*ones(1,N+1)];
    otherwise
        error('undefined reference')
end

% convert reference to Pi-space
xiref = C * Mxinv * C' * xiref;

%% simulation

% linearization data from the previous step
oldXis = repmat(xi0,1,N+1);
oldUs = repmat(u0,1,M);
xihats = nan(nx,N+1);

% Jacobian matrices
As = cell(1,N);
Bs = cell(1,N);

% log data
Xis = nan(nx,nsim+1);
Us = nan(nu,nsim);
Xis(:,1) = xi0;

% figure; hold on  % for live plotting
Xi = xi0;
U = u0;
ax = 0;
ay = 0;

tic
for i = 1 : nsim
    % linearize around the predicted states and inputs
    xihat = Xi;
    xihats(:,1) = Mxinv * Xi;
    for j = 1:N
        if j > M
            [~,xihat] = ode45(@OdeTwoTrackModel,[0,Ts],[xihat;oldUs(:,M);deltaf(i)]);
            [A,B] = TwoTrackModelJacobian(oldXis(:,j),oldUs(:,M),deltaf(i));
        else
            [~,xihat] = ode45(@OdeTwoTrackModel,[0,Ts],[xihat;oldUs(:,j);deltaf(i)]);
            [A,B] = TwoTrackModelJacobian(oldXis(:,j),oldUs(:,j),deltaf(i));
        end
        xihat = xihat(end,1:nx)';  % nx*(N+1)
        % transform to Pi-space
        xihats(:,j+1) = Mxinv * xihat;
        
        % formulate in Pi-space and discretize
        Ap = that * Mxinv * A * Mxi;
        Bp = that * Mxinv * B * Mu;
        [Ad,Bd] = c2d(Ap,Bp,Tsp);
        As{j} = Ad;
        Bs{j} = Bd;       
    end
    
    % convert current state and past inputs to Pi-space
    Xi = Mxinv * Xi;
    oldUs = Muinv * oldUs;
    
    % solve the optimization problem
    inputs = {Xi,oldUs,xiref(:,i:i+N),[As{:}],[Bs{:}],xihats,deltaf(i)};
    [solutions, errorcode, errortext, ~, ~, diagnostics] = controller(inputs);
    
    % check the optimization results
    if errorcode ~= 0
        error(errortext{1});
    end
    
    % "live" plot
    predictions = solutions{2};
%     cla;
%     stairs(i:i+N,Mxi(2,2)*xiref(2,i:i+N),'k')
%     stairs(i:i+N,Mxi(2,2)*[Xi(2),predictions(2,2:end)],'b')
%     stairs(1:i,Xis(2,1:i),'g');title('Prediction+reference+state')
%     xlabel('k')
%     ylabel('vx [m/s')
%     legend('reference','prediction','trajectory')
%     pause(0.01)
    
    % apply the inputs to the plant
    U = solutions{1}(:,1);  % take the first optimal input
    U = Mu * U  % convert back from Pi-space
    Xi = Mxi * Xi;
    [~,Xi] = ode45(@OdeSimulationModel,[0,Ts],[Xi;U;deltaf(i);ax;ay]);
    Xi = Xi(end,1:nx)'
    i
    
    % calculate the acceleration
    f = OdeSimulationModel(0,[Xi;U;deltaf(i);ax;ay]);
    ax = f(1); ay = f(2);
    
    % save optimization results for the next iteration
    oldXis = [Xi, Mxi * predictions(:,2:end)];
    oldUs = Mu * solutions{1};
    
    % save data for plotting
    Xis(:,i+1) = Xi;
    Us(:,i) = U;
end
toc

%% plot the results
t = 0:Ts:nsim*Ts;
xiref = Mxi(1:3,1:3) * xiref(:,1:nsim+1);

% states and reference
figure
ylabels = {'$v_x$ [m/s]','$v_y$ [m/s]','$\dot{\theta}$ [rad/s]',...
           '\omega FL [rad/s]','\omega FR [rad/s]','\omega RL [rad/s]','\omega RR [rad/s]'};
for i = 1:3
    subplot(3,1,i)
    plot(t,Xis(i,:))
    hold on
    plot(t,xiref(i,:))
    xlabel('Time [s]')
    ylabel(ylabels{i})
end
subplot(3,1,1)
title('States and reference')

% inputs
figure
stairs(t,[u0 Us]')
legend('T_{fl}','T_{fr}','T_{rl}','T_{rr}')
title('Input torque')
xlabel('t [s]')
ylabel('Input torque [Nm]')

% input rate
figure
stairs(t(2:end),diff([u0 Us],1,2)')
legend('\Delta T_{fl}','\Delta T_{fr}','\Delta T_{rl}','\Delta T_{rr}')
title('Input torque rate')
xlabel('t [s]')
ylabel('Input torque rate [Nm/s]')

% plot the trajectory
dx = Xis(1,:);
dy = Xis(2,:);
psi = Ts*cumtrapz(Xis(3,:));
dX = dx.*cos(psi) - dy.*sin(psi);
dY = dx.*sin(psi) + dy.*cos(psi);
X = Ts*cumtrapz(dX);
Y = Ts*cumtrapz(dY);

X1 =  lf*cos(psi) - w/2*sin(psi);
X2 =  lf*cos(psi) + w/2*sin(psi);
X3 = -lr*cos(psi) + w/2*sin(psi);
X4 = -lr*cos(psi) - w/2*sin(psi);

Y1 =  w/2*cos(psi) + lf*sin(psi);
Y2 = -w/2*cos(psi) + lf*sin(psi);
Y3 = -w/2*cos(psi) - lf*sin(psi);
Y4 =  w/2*cos(psi) - lf*sin(psi);

figure
plot(X,Y,'b--')
axis equal;
grid on;
hold on;
title('Path')
xlabel('X [m]')
ylabel('Y [m]')

for i = 1:1/(Ts):length(psi)
    patch(X(i)+[X1(i) X2(i) X3(i) X4(i)], Y(i)+[Y1(i) Y2(i) Y3(i) Y4(i)],'red','facecolor','none','edgecolor','red');
end

%% old code archive

%     slipConstraints = [    
%         params.rw*xi{k}(5) - (slipMax + 1)*(cos(deltaf)*(xi{k}(2) - xi{k}(4)*params.w) + sin(deltaf)*(xi{k}(1) + params.lf*xi{k}(4))) <= e{k}(1),...
%         params.rw*xi{k}(6) - (slipMax + 1)*(cos(deltaf)*(xi{k}(2) + xi{k}(4)*params.w) + sin(deltaf)*(xi{k}(1) + params.lf*xi{k}(4))) <= e{k}(2),...
%         params.rw*xi{k}(7) - (xi{k}(2) - xi{k}(4)*params.w)*(slipMax + 1) <= e{k}(3),...
%         params.rw*xi{k}(8) - (xi{k}(2) + xi{k}(4)*params.w)*(slipMax + 1) <= e{k}(4),...
%         (slipMin + 1)*(cos(deltaf)*(xi{k}(2) - xi{k}(4)*params.w) + sin(deltaf)*(xi{k}(1) + params.lf*xi{k}(4))) - params.rw*xi{k}(5) <= e{k}(5),...
%         (slipMin + 1)*(cos(deltaf)*(xi{k}(2) + xi{k}(4)*params.w) + sin(deltaf)*(xi{k}(1) + params.lf*xi{k}(4))) - params.rw*xi{k}(6) <= e{k}(6),...
%         (slipMin + 1)*(xi{k}(2) - xi{k}(4)*params.w) - params.rw*xi{k}(7) <= e{k}(7),...
%         (slipMin + 1)*(xi{k}(2) + xi{k}(4)*params.w) - params.rw*xi{k}(8) <= e{k}(8)]

%     slipAngleConstraints = [        
%         cos(deltaf)*(xi{k}(1) + params.lf*xi{k}(4)) - sin(deltaf)*(xi{k}(2) - xi{k}(4)*params.w) - slipAngleMax*(cos(deltaf)*(xi{k}(2) - xi{k}(4)*params.w) + sin(deltaf)*(xi{k}(1) + params.lf*xi{k}(4))) <= e{k}(9),...
%         cos(deltaf)*(xi{k}(1) + params.lf*xi{k}(4)) - sin(deltaf)*(xi{k}(2) + xi{k}(4)*params.w) - slipAngleMax*(cos(deltaf)*(xi{k}(2) + xi{k}(4)*params.w) + sin(deltaf)*(xi{k}(1) + params.lf*xi{k}(4))) <= e{k}(10),...
%         xi{k}(1) - params.lr*xi{k}(4) - slipAngleMax*(xi{k}(2) - xi{k}(4)*params.w) <= e{k}(11),...
%         xi{k}(1) - params.lr*xi{k}(4) - slipAngleMax*(xi{k}(2) + xi{k}(4)*params.w) <= e{k}(12),...
%         slipAngleMin*(cos(deltaf)*(xi{k}(2) - xi{k}(4)*params.w) + sin(deltaf)*(xi{k}(1) + params.lf*xi{k}(4))) - cos(deltaf)*(xi{k}(1) + params.lf*xi{k}(4)) - sin(deltaf)*(xi{k}(2) - xi{k}(4)*params.w) <= e{k}(13),...
%         slipAngleMin*(cos(deltaf)*(xi{k}(2) + xi{k}(4)*params.w) + sin(deltaf)*(xi{k}(1) + params.lf*xi{k}(4))) - cos(deltaf)*(xi{k}(1) + params.lf*xi{k}(4)) - sin(deltaf)*(xi{k}(2) + xi{k}(4)*params.w) <= e{k}(14),...
%         slipAngleMin*(xi{k}(2) - xi{k}(4)*params.w) - xi{k}(1) - params.lr*xi{k}(4) <= e{k}(15),...
%         slipAngleMin*(xi{k}(2) + xi{k}(4)*params.w) - xi{k}(1) - params.lr*xi{k}(4) <= e{k}(16)];

% plot(t,Xis')
% hold on
% plot(t,(Mxi(1:4,1:4)*xiref(:,1:nsim+1))')
% legend('v_y','v_x','\theta','der(\theta)',...
%        '\omega_{fl}','\omega_{fr}','\omega_{rl}','\omega_{rr}',...
%        'v_{y,ref}','v_{x,ref}','\theta_{ref}','der(\theta)_{ref}')
% xlabel('t [s]')
% ylabel('States and reference')