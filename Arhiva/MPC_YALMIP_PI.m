%% Two-track model LTV MPC
%% problem setup
clear;yalmip('clear');close all;clc
load params.mat  % vehicle parameters
nx = 8;  % no. of states
nu = 4;  % no. of inputs
C = [eye(4),zeros(4)];
ny = size(C,1);  % no. of tracked states
Ts = 0.05;  % sampling time
N = 10;  % horizon length
M = 3;  % control horizon (blocking)
nsim = 102;  % simulation length

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

% YALMIP data
u = sdpvar(nu*ones(1,N),ones(1,N));  % inputs
xi = sdpvar(nx*ones(1,N+1),ones(1,N+1));  % states
e = sdpvar(16*ones(1,N+1),ones(1,N+1));  % soft constraint variables
pastu = sdpvar(4,1);  % for input rate limits
ops = sdpsettings('verbose',1,...
                  'solver','osqp',...
                  'osqp.eps_abs',1e-6,...
                  'osqp.eps_rel',1e-6,...
                  'debug',1,...
                  'cachesolvers',1);  % print output
ops.savesolveroutput = 1;
ops.savesolverinput = 1;
% ops.usex0 = 1;
% ops = sdpsettings('solver','quadprog');

% QP weights
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

% initial state and input
vx0 = 10;
xi0 = [0; vx0; 0; 0; 250/79*vx0*ones(4,1)];
u0 = 3.2590*ones(4,1);
u0 = zeros(4,1);

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
                 [vx0*ones(1,nsim/3),vx0+(vx1-vx0)/nsim*3*(1:nsim/3),vx1*ones(1,nsim/3+N+1)];
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

% convert reference to Pi-space
xirefmat = C * Mxinv * C' * xiref;
xiref = num2cell(xirefmat,1);  % split columnwise into cell array

% linearization data from the previous step
oldXis = repmat(xi0,1,N+1);
oldUs = repmat(u0,1,M);
xihats = cell(1,N+1);

% Jacobian matrices
As = cell(1,N);
Bs = cell(1,N);

% log data
Xis = nan(nx,nsim+1);
Us = nan(nu,nsim);
Xis(:,1) = xi0;

%% simulation
figure; hold on
Xi = xi0;
U = u0;
ay = 0;
ax = 0;

tic
for i = 1 : nsim    
    % linearize around the predicted states and inputs
    xihat = Xi;
    xihats{1} = Mxinv * Xi;
    for j = 1:N
        if j > M
            [~,xihat] = ode45(@odeTwoTrackModel,[0,Ts],[xihat;oldUs(:,M);deltaf(i)]);
            [A,B] = twoTrackModelJacobian(oldXis(:,j),oldUs(:,M),deltaf(i));
        else
            [~,xihat] = ode45(@odeTwoTrackModel,[0,Ts],[xihat;oldUs(:,j);deltaf(i)]);
            [A,B] = twoTrackModelJacobian(oldXis(:,j),oldUs(:,j),deltaf(i));
        end
        xihat = xihat(end,1:8)';  % nx*(N+1)
        % transform to Pi-space
        xihats{j+1} = Mxinv * xihat;
        
        % formulate in Pi-space and discretize
        Ap = that * Mxinv * A * Mxi;
        Bp = that * Mxinv * B * Mu;
        [Ad,Bd] = c2d(Ap,Bp,Tsp);
        As{j} = Ad;
        Bs{j} = Bd;
    end
    
    % convert current state and past inputs to Pi-space
    Xi = Mxinv * Xi;
    U = Muinv * U;
    oldUs = Muinv * oldUs;
    
    % formulate the optimization problem
    deltau = diff([pastu u{:}],1,2);
    constraints = [(xi{1} == Xi):'Initial state'];  % initial state constraint
    constraints = [constraints, (pastu == U):'Past input'];
    constraints = [constraints, (-torqueMax <= [u{:}] <= torqueMax):'Input limits'];
    constraints = [constraints, (-torqueRateMax <= deltau <= torqueRateMax):'Input rate limits'];
    constraints = [constraints, ([e{:}] >= 0):'Positive slack'];
    objective = 0;
    for k = 1:N+1
        % system dynamics and the objective function
        if k == N+1
            objective = objective + (C*xi{k}-xiref{i+k-1})'*QN*(C*xi{k}-xiref{i+k-1});
        else
            if k > M
                objective = objective + u{M}'*S*u{M};
                constraints = [constraints, (xi{k+1} == As{k}*xi{k} + Bs{k}*u{M} + xihats{k+1} - As{k}*xihats{k} - Bs{k}*oldUs(:,M)):['Dynamics ' num2str(k)]];
            else
                objective = objective + u{k}'*S*u{k};
                constraints = [constraints, (xi{k+1} == As{k}*xi{k} + Bs{k}*u{k} + xihats{k+1} - As{k}*xihats{k} - Bs{k}*oldUs(:,k)):['Dynamics ' num2str(k)]];
            end
            objective = objective + (C*xi{k}-xiref{i+k-1})'*Q*(C*xi{k}-xiref{i+k-1});
            if k < M
                objective = objective + deltau(:,k)' * R * deltau(:,k);
            end
        end
        
        % state constraints
        % tire slip constraints
        slipConstraints = [    
            params.rw*xi{k}(5) - (slipMax + 1)*(cos(deltaf(i))*(xi{k}(2) - xi{k}(4)*params.w) + sin(deltaf(i))*(xi{k}(1) + params.lf*xi{k}(4))) <= e{k}(1),...
            params.rw*xi{k}(6) - (slipMax + 1)*(cos(deltaf(i))*(xi{k}(2) + xi{k}(4)*params.w) + sin(deltaf(i))*(xi{k}(1) + params.lf*xi{k}(4))) <= e{k}(2),...
            params.rw*xi{k}(7) - (xi{k}(2) - xi{k}(4)*params.w)*(slipMax + 1) <= e{k}(3),...
            params.rw*xi{k}(8) - (xi{k}(2) - xi{k}(4)*params.w)*(slipMax + 1) <= e{k}(4),...
            -(slipMax - 1)*(cos(deltaf(i))*(xi{k}(2) - xi{k}(4)*params.w) + sin(deltaf(i))*(xi{k}(1) + params.lf*xi{k}(4))) - params.rw*xi{k}(5) <= e{k}(5),...
            -(slipMax - 1)*(cos(deltaf(i))*(xi{k}(2) + xi{k}(4)*params.w) + sin(deltaf(i))*(xi{k}(1) + params.lf*xi{k}(4))) - params.rw*xi{k}(6) <= e{k}(6),...
            -(slipMax - 1)*(xi{k}(2) - xi{k}(4)*params.w) - params.rw*xi{k}(7) <= e{k}(7),...
            -(slipMax - 1)*(xi{k}(2) - xi{k}(4)*params.w) - params.rw*xi{k}(8) <= e{k}(8)];

        % tire slip angle constraints
        slipAngleConstraints = [        
            cos(deltaf(i))*(xi{k}(1) + params.lf*xi{k}(4)) - sin(deltaf(i))*(xi{k}(2) - xi{k}(4)*params.w) - slipAngleMax*(cos(deltaf(i))*(xi{k}(2) - xi{k}(4)*params.w) + sin(deltaf(i))*(xi{k}(1) + params.lf*xi{k}(4))) <= e{k}(9),...
            cos(deltaf(i))*(xi{k}(1) + params.lf*xi{k}(4)) - sin(deltaf(i))*(xi{k}(2) + xi{k}(4)*params.w) - slipAngleMax*(cos(deltaf(i))*(xi{k}(2) + xi{k}(4)*params.w) + sin(deltaf(i))*(xi{k}(1) + params.lf*xi{k}(4))) <= e{k}(10),...
            xi{k}(1) - params.lr*xi{k}(4) - slipAngleMax*(xi{k}(2) - xi{k}(4)*params.w) <= e{k}(11),...
            xi{k}(1) - params.lr*xi{k}(4) - slipAngleMax*(xi{k}(2) + xi{k}(4)*params.w) <= e{k}(12),...
            -slipAngleMax*(cos(deltaf(i))*(xi{k}(2) - xi{k}(4)*params.w) + sin(deltaf(i))*(xi{k}(1) + params.lf*xi{k}(4))) - cos(deltaf(i))*(xi{k}(1) + params.lf*xi{k}(4)) - sin(deltaf(i))*(xi{k}(2) - xi{k}(4)*params.w) <= e{k}(13),...
            -slipAngleMax*(cos(deltaf(i))*(xi{k}(2) + xi{k}(4)*params.w) + sin(deltaf(i))*(xi{k}(1) + params.lf*xi{k}(4))) - cos(deltaf(i))*(xi{k}(1) + params.lf*xi{k}(4)) - sin(deltaf(i))*(xi{k}(2) + xi{k}(4)*params.w) <= e{k}(14),...
            -slipAngleMax*(xi{k}(2) - xi{k}(4)*params.w) - xi{k}(1) - params.lr*xi{k}(4) <= e{k}(15),...
            -slipAngleMax*(xi{k}(2) + xi{k}(4)*params.w) - xi{k}(1) - params.lr*xi{k}(4) <= e{k}(16)];

        constraints = [constraints, slipConstraints:['Slip ' num2str(k)],...
                                    slipAngleConstraints:['Slip angle ' num2str(k)]];

        objective = objective + p * e{k}'*e{k};  % penalize slack variables     
    end
    
    % solve the optimization problem
    diagnostics = optimize(constraints,objective,ops);
    
    % check the optimization results
    if diagnostics.problem ~= 0
        error(diagnostics.info);
    end
    
    % "live" plot
    predictions = value([xi{:}]);
    cla;
    stairs(i:i+N,Mxi(2,2)*xirefmat(2,i:i+N),'k')
    stairs(i:i+N,Mxi(2,2)*predictions(2,:),'b')
    stairs(1:i,Xis(2,1:i),'g');title('Prediction+reference+state')
    xlabel('k')
    ylabel('vx [m/s')
    legend('reference','prediction','trajectory')
    pause(0.01)
    
    % Apply the nu inputs to the plant
    U = value(u{1});  % take the first optimal input
    U = Mu * U  % convert back from Pi-space
    Xi = Mxi * Xi;
    
    [~,Xi] = ode45(@odeSimulationModel,[0,Ts],[Xi;U;deltaf(i);ay;ax]);
    Xi = Xi(end,1:8)'
    i
    
    % calculate the acceleration
    f = odeSimulationModel(0,[Xi;U;deltaf(i);ay;ax]);
    ay = f(1); ax = f(2);
    
    % save optimization results for the next iteration
    oldXis = [Xi, Mxi * predictions(:,2:end)];
    oldUs = Mu * value([u{:}]);
    
    % save data for plotting
    Xis(:,i+1) = Xi;
    Us(:,i) = U;
end
toc

%% plot the results
t = 0:Ts:nsim*Ts;
xiref = Mxi(1:4,1:4) * xiref(:,1:nsim+1);

% states and reference
figure
ylabels = {'vy [m/s]','vx [m/s]','\theta [rad]','d\theta [rad/s]',...
           '\omega FL [rad/s]','\omega FR [rad/s]','\omega RL [rad/s]','\omega RR [rad/s]'};
for i = 1:4
    subplot(4,1,i)
    plot(t,Xis(i,:))
    hold on
    plot(t,xiref(i,:))
    xlabel('Time [s]')
    ylabel(ylabels{i})
end
subplot(4,1,1)
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
dy = Xis(1,:);
dx = Xis(2,:);
psi = Xis(3,:);
dY = dx.*sin(psi) + dy.*cos(psi);
dX = dx.*cos(psi) - dy.*sin(psi);
Y = Ts*cumtrapz(dY);
X = Ts*cumtrapz(dX);

lf = params.lf;
lr = params.lr;
w  = params.w;

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

%%
%         [~,xihat] = ode45(@odeTwoTrackModel,[0,Ts],[xihat;oldUs(:,j);deltaf(i+j-1)]);
%         xihat = xihat(end,1:8)';
%         % transform to Pi-space
%         xihats(:,j+1) = Mxinv * xihat;
%         
%         [A,B] = twoTrackModelJacobian(oldXis(:,j),oldUs(:,j),deltaf(i+j-1));