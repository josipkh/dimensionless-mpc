%% Two-track model LTV MPC

clear;clc;close all

generate_code = false;
generated_code_folder = './Generated code';

% vehicle parameters
load Params.mat

% system dimensions
nx = params.nx;             % no. of states
nu = params.nu;             % no. of inputs
neps = 16;                  % no. of slack variables
C = [eye(4),zeros(4)];      % first four states are available
ny = size(C,1);             % no. of tracked states

% simulation setup
Ts = 0.05;                  % sampling time [s]
Tsim = 5;                   % simulation time [s]
nsim = Tsim / Ts;           % no. of simulation steps

% MPC parameters
N = 10;                                 % prediction horizon
M = 3;                                  % control horizon (blocking)
Q = diag([10 100 1e4 2e4 0 0 0 0]);     % tracked state weights
QN = Q;                                 % terminal cost weights
S = 1e-3*eye(nu);                       % input weights
R = 1e-2*eye(nu);                       % input rate weights
p = 100;                                % slack variable weights

enablePi = 1;       % use Pi-groups?
constraints = [1;   % dynamics
               0;   % state (slip, slip angle) -> soft constraints
               1;   % input
               0];  % input rate

% Pi-groups
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
% Cp = C * Mxi;

% transform MPC weights to Pi-space
Q = Mxi * Q * Mxi;
S = Mu * S * Mu;
R = that^-2 * Mu * R * Mu;

% initialize optimization problem matrices
nopt = nx*(N+1) + nu*M + nu*M + neps*(N+1);
A = [];
l = [];
u = [];

%% Objective function
% Cast MPC problem to a QP: x = (x(0),x(1),...,x(N),u(0),...,u(M-1))
% optimization vector size: nx*(N+1) + nu*M + nu*M + neps*(N+1)

% quadratic objective
P = blkdiag( kron(speye(N), Q), QN, kron(speye(M), S), ...
             kron(speye(M), R), kron(speye(N+1), p*eye(neps)));
         
% linear objective
xiref = rand(nx,N);
q = [reshape(-Q*xiref, N*nx, 1); -QN*xiref(:,end); 
     zeros(M*nu, 1); zeros(M*nu, 1); zeros((N+1)*neps, 1)];
 
%% system dynamics
As = mat2cell(rand(nx,nx*N), nx, nx*ones(1,N));
Bs = mat2cell(rand(nx,nu*N), nx, nu*ones(1,N));

% Ax = [-1   0   0   0
%       A1  -1   0   0
%        0  A2  -1   0
%        0   0  AN  -1]
Ax = kron(speye(N+1), -speye(nx)) + ...
     [sparse(nx, nx*(N+1)); blkdiag(As{:}), zeros(nx*N, nx)];
 
% Bu = [0   0 
%       B1  0 
%       0  BM
%       0  BN] (blocking)
blockingBlock = cell2mat(Bs(M+1:end)');  % check this block to remove warnings
if isempty(blockingBlock); blockingBlock = zeros(nx*(N-M),nu); end
Bu = [sparse(nx, nu*M); 
      blkdiag(Bs{1:M});
      zeros(nx*(N-M),nu*(M-1)), blockingBlock];

ADynamics = [Ax, Bu, zeros(size(Ax,1), nu*M + neps*(N+1))];

Xi = rand(nx,1);  % initial state
d = rand(nx,N);   % linearization offset
lDynamics = [-Xi; -reshape(d, N*nx, 1)];
uDynamics = lDynamics;

if constraints(1)
    A = [A; ADynamics];
    l = [l; lDynamics];
    u = [u; uDynamics];
end

%% state constraints (soft)

% slip and slip angle limits (for stability)
slipMax = params.slipMax;  % conversion not needed (dimensionless)
slipMin = params.slipMin;
% slipMax = 0.3; slipMin = -0.3;
alphaMax = params.alphaMax;  % conversion not needed (angle)
alphaMin = params.alphaMin;
% alphaMax = 0.1745; alphaMin = -0.1745;
w = params.w; rw = params.rw; lf = params.lf; lr = params.lr;
if enablePi; le = lf+lr; w = w/le; rw = rw/le; lf = lf/le; lr = lr/le; end
deltaf = deg2rad(5);

slipFLUpperBound = [
    (-sin(deltaf)*(slipMax + 1)), (-cos(deltaf)*(slipMax + 1)), 0,...
    (w*cos(deltaf) - lf*sin(deltaf))*(slipMax + 1), rw, 0, 0, 0];
slipFRUpperBound = [
    (-sin(deltaf)*(slipMax + 1)), (-cos(deltaf)*(slipMax + 1)), 0,...
    (-(w*cos(deltaf) + lf*sin(deltaf))*(slipMax + 1)), 0, rw, 0, 0];
slipRLUpperBound = [0, (-slipMax - 1), 0, w*(slipMax + 1), 0, 0, rw, 0];
slipRRUpperBound = [0, (-slipMax - 1), 0, (-w*(slipMax + 1)), 0, 0, 0, rw];

slipFLLowerBound = [
    sin(deltaf)*(slipMin + 1), cos(deltaf)*(slipMin + 1), 0,...
    (-(w*cos(deltaf) - lf*sin(deltaf))*(slipMin + 1)), -rw, 0, 0, 0];
slipFRLowerBound = [
    sin(deltaf)*(slipMin + 1), cos(deltaf)*(slipMin + 1), 0,...
    (w*cos(deltaf) + lf*sin(deltaf))*(slipMin + 1), 0, -rw, 0, 0];
slipRLLowerBound = [0, (slipMin + 1), 0, -w*(slipMin + 1), 0, 0, -rw, 0];
slipRRLowerBound = [0, (slipMin + 1), 0, w*(slipMin + 1), 0, 0, 0, -rw];

alphaFLUpperBound = [
    (cos(deltaf) - alphaMax*sin(deltaf)), (- sin(deltaf) - alphaMax*cos(deltaf)), 0,...
    (alphaMax*(w*cos(deltaf) - lf*sin(deltaf)) + lf*cos(deltaf) + w*sin(deltaf)), zeros(1,4)];
alphaFRUpperBound = [
    (cos(deltaf) - alphaMax*sin(deltaf)), (- sin(deltaf) - alphaMax*cos(deltaf)), 0,...
    (lf*cos(deltaf) - alphaMax*(w*cos(deltaf) + lf*sin(deltaf)) - w*sin(deltaf)), zeros(1,4)];
alphaRLUpperBound = [1, (-alphaMax), 0, (alphaMax*w - lr), zeros(1,4)];
alphaRRUpperBound = [1, (-alphaMax), 0, (- lr - alphaMax*w), zeros(1,4)];

alphaFLLowerBound = [
    (alphaMin*sin(deltaf) - cos(deltaf)), (sin(deltaf) + alphaMin*cos(deltaf)), 0,...
    (- alphaMin*(w*cos(deltaf) - lf*sin(deltaf)) - lf*cos(deltaf) - w*sin(deltaf)), zeros(1,4)];
alphaFRLowerBound = [
    (alphaMin*sin(deltaf) - cos(deltaf)), (sin(deltaf) + alphaMin*cos(deltaf)), 0,...
    (alphaMin*(w*cos(deltaf) + lf*sin(deltaf)) - lf*cos(deltaf) + w*sin(deltaf)), zeros(1,4)];
alphaRLLowerBound = [-1, alphaMin, 0, lr - alphaMin*w, zeros(1,4)];
alphaRRLowerBound = [-1, alphaMin, 0, lr + alphaMin*w, zeros(1,4)];

stateConstraintBlock = [
    slipFLUpperBound; slipFRUpperBound; slipRLUpperBound; slipRRUpperBound;...
    slipFLLowerBound; slipFRLowerBound; slipRLLowerBound; slipRRLowerBound;...
    alphaFLUpperBound; alphaFRUpperBound; alphaRLUpperBound; alphaRRUpperBound;...
    alphaFLLowerBound; alphaFRLowerBound; alphaRLLowerBound; alphaRRLowerBound];

stateConstraintBlock = stateConstraintBlock * Mxi;

AStateConstraint = [kron(speye(N+1), stateConstraintBlock), zeros(neps*(N+1),2*nu*M), -speye(neps*(N+1))];
lStateConstraint = -inf(size(AStateConstraint,1),1);
uStateConstraint = zeros(size(AStateConstraint,1),1);

if constraints(2)
    A = [A; AStateConstraint];
    l = [l; lStateConstraint];
    u = [u; uStateConstraint];
end

%% actuator limitations

% input constraints
torqueMax = params.Tmax;  % N
torqueMax = Muinv(1) * torqueMax;  % convert to Pi-space
torqueMin = Muinv(1) * params.Tmin;
% torqueMax = 1000; torqueMin = -1000;

AInputConstraint = [zeros(nu*M,nx*(N+1)), speye(nu*M), zeros(nu*M,nu*M+neps*(N+1))];
lInputConstraint = torqueMin * ones(nu*M,1);
uInputConstraint = torqueMax * ones(nu*M,1);

if constraints(3)
    A = [A; AInputConstraint];
    l = [l; lInputConstraint];
    u = [u; uInputConstraint];
end   

% input rate constraints
torqueRateMax = params.dTmax;  % Nm/s
torqueRateMax = Muinv(1) * that * torqueRateMax;  % convert to Pi-space
torqueRateMax = torqueRateMax * Tsp;  % max change in one time step
torqueRateMin = Muinv(1) * that * params.dTmin * Tsp;
% torqueRateMax = torqueRateMax / that / Tsp * Ts;
% torqueRateMin = torqueRateMin / that / Tsp * Ts;

AInputRateConstraint = [zeros(nu*M,nx*(N+1)+nu*M), speye(nu*M), zeros(nu*M,neps*(N+1))];
lInputRateConstraint = torqueRateMin * ones(nu*M,1);
uInputRateConstraint = torqueRateMax * ones(nu*M,1);

if constraints(4)
    A = [A; AInputRateConstraint];
    l = [l; lInputRateConstraint];
    u = [u; uInputRateConstraint];
end

%% other (in)equalities

% deltaU(k) = u(k) - u(k-1);  u(-1) = previous input
% -u0 <= 0 ... 0 -u1 0 ... 0 du1 0 ... 0 <= -u0
%  0  <= 0 ... 0  u(k-1) -u(k) 0 ... 0 du(k) ... 0 <= 0
AdeltaU = [zeros(nu*M,nx*(N+1)), -speye(nu*M) + diag(ones(nu*(M-1),1),-nu), speye(nu*M), zeros(nu*M,neps*(N+1))];
ldeltaU = [-rand(nu,1); zeros(nu*(M-1),1)];
udeltaU = ldeltaU;

if constraints(4)
    A = [A; AdeltaU];
    l = [l; ldeltaU];
    u = [u; udeltaU];
end

% slack variables are positive
ASlack = [zeros(neps*(N+1), nx*(N+1)+2*nu*M), speye(neps*(N+1))];
lSlack = 0*ones(neps*(N+1),1);
uSlack = inf*ones(neps*(N+1),1);

if constraints(2)
    A = [A; ASlack];
    l = [l; lSlack];
    u = [u; uSlack];
end

%% Create the controller

prob = osqp;
prob.setup(P, q, A, l, u,...
           'warm_start', true,...
           'eps_abs',1e-6,...
           'eps_rel',1e-6,...
           'verbose',0);
       
% save the indices for later updating (problem when zeros appear)
Ax_idx = find(A);

% generate the code if needed
if generate_code
    prob.codegen(generated_code_folder,'parameters','matrices')
end

%% initial state and reference

% initial state and input
vx0 = 10;   % [m/s] initial long. speed
xi0 = [0; vx0; 0; 0; 250/79*vx0*ones(4,1)];
u0 = 3.2590*ones(4,1);

% reference generation
reference = 2;
switch reference
    case 1
        % constant forward driving
        u0 = zeros(4,1);
        deltaf = zeros(1,nsim);
        xiref = repmat(C*xi0,1,nsim+N+1);
    case 2
        % change speed from vx0 to vx1
        u0 = zeros(4,1);
        deltaf = zeros(1,nsim);
        vx1 = 20;
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

% add reference for untracked states to respect the OSQP format
xiref = [xiref; zeros(nx-ny,size(xiref,2))];
% convert reference to Pi-space
xiref = Mxinv * xiref;

%% Simulate in closed loop

% linearization data from the previous step
oldXis = repmat(xi0,1,N+1);
oldUs = repmat(u0,1,M);
xihats = nan(nx,N+1);
d = nan(nx,N);

% Jacobian matrices
As = cell(1,N);
Bs = cell(1,N);

% log data
Xis = nan(nx,nsim+1);
Us = nan(nu,nsim);
Xis(:,1) = xi0;

% set the initial state
Xi = xi0;
U = u0;
ay = 0;
ax = 0;
% figure; hold on

tic
for i = 1 : nsim
    % linearize around the predicted states and inputs
    xihat = Xi;  % start at the current system state (normal space)
    xihats(:,1) = Mxinv * Xi;  % save predictions in Pi-space
    for j = 1:N
        if j > M
            % integrate the model to get the predictions
            [~,xihat] = ode45(@OdeTwoTrackModel,[0,Ts],[xihat;oldUs(:,M);deltaf(i)]);
            [A,B] = TwoTrackModelJacobian(oldXis(:,j),oldUs(:,M),deltaf(i));
        else
            [~,xihat] = ode45(@OdeTwoTrackModel,[0,Ts],[xihat;oldUs(:,j);deltaf(i)]);
            [A,B] = TwoTrackModelJacobian(oldXis(:,j),oldUs(:,j),deltaf(i));
        end
        xihat = xihat(end,1:8)';  % take the value at Ts
        xihats(:,j+1) = Mxinv * xihat;  % transform to Pi-space
        
        % formulate in Pi-space and discretize
        Ap = that * Mxinv * A * Mxi;
        Bp = that * Mxinv * B * Mu;
        [Ad,Bd] = c2d(Ap,Bp,Tsp);
        As{j} = Ad;
        Bs{j} = Bd; 
        
        % save the linearization offset
        if j > M
            d(:,j) = xihats(:,j+1) - Ad * xihats(:,j) - Bd * Muinv * oldUs(:,M);  % convert the past input
        else
            d(:,j) = xihats(:,j+1) - Ad * xihats(:,j) - Bd * Muinv * oldUs(:,j);
        end
    end
    
    % convert current state to Pi-space
    Xi = Mxinv * Xi;
    
    % update the linear cost
    q = [reshape(-Q*xiref(:,i:i+N-1), N*nx, 1); -QN*xiref(:,i+N); 
         zeros(M*nu, 1); zeros(M*nu, 1); zeros((N+1)*neps, 1)];
    
    % update the system dynamics constraints
    Ax = kron(speye(N+1), -speye(nx)) + [sparse(nx, nx*(N+1)); blkdiag(As{:}), zeros(nx*N, nx)];
    if M < N; blockingBlock = cell2mat(Bs(M+1:end)'); end  % otherwise zeros from setup
    Bu = [sparse(nx, nu*M); blkdiag(Bs{1:M}); zeros(nx*(N-M),nu*(M-1)), blockingBlock];
    ADynamics = [Ax, Bu, zeros(size(Ax,1), nu*M + neps*(N+1))];
    lDynamics = [-Xi; -reshape(d, N*nx, 1)];
    uDynamics = lDynamics;
    
    % update the input rate constraints    
    ldeltaU = [-Muinv * U; zeros(nu*(M-1),1)];
    udeltaU = ldeltaU;
    
    % update the state constraints
    slipFLUpperBound = [
        (-sin(deltaf(i))*(slipMax + 1)), (-cos(deltaf(i))*(slipMax + 1)), 0,...
        (w*cos(deltaf(i)) - lf*sin(deltaf(i)))*(slipMax + 1), rw, 0, 0, 0];
    slipFRUpperBound = [
        (-sin(deltaf(i))*(slipMax + 1)), (-cos(deltaf(i))*(slipMax + 1)), 0,...
        (-(w*cos(deltaf(i)) + lf*sin(deltaf(i)))*(slipMax + 1)), 0, rw, 0, 0];
    slipFLLowerBound = [
        sin(deltaf(i))*(slipMin + 1), cos(deltaf(i))*(slipMin + 1), 0,...
        (-(w*cos(deltaf(i)) - lf*sin(deltaf(i)))*(slipMin + 1)), -rw, 0, 0, 0];
    slipFRLowerBound = [
        sin(deltaf(i))*(slipMin + 1), cos(deltaf(i))*(slipMin + 1), 0,...
        (w*cos(deltaf(i)) + lf*sin(deltaf(i)))*(slipMin + 1), 0, -rw, 0, 0];
    alphaFLUpperBound = [
        (cos(deltaf(i)) - alphaMax*sin(deltaf(i))), (- sin(deltaf(i)) - alphaMax*cos(deltaf(i))), 0,...
        (alphaMax*(w*cos(deltaf(i)) - lf*sin(deltaf(i))) + lf*cos(deltaf(i)) + w*sin(deltaf(i))), zeros(1,4)];
    alphaFRUpperBound = [
        (cos(deltaf(i)) - alphaMax*sin(deltaf(i))), (- sin(deltaf(i)) - alphaMax*cos(deltaf(i))), 0,...
        (lf*cos(deltaf(i)) - alphaMax*(w*cos(deltaf(i)) + lf*sin(deltaf(i))) - w*sin(deltaf(i))), zeros(1,4)];
    alphaFLLowerBound = [
        (alphaMin*sin(deltaf(i)) - cos(deltaf(i))), (sin(deltaf(i)) + alphaMin*cos(deltaf(i))), 0,...
        (- alphaMin*(w*cos(deltaf(i)) - lf*sin(deltaf(i))) - lf*cos(deltaf(i)) - w*sin(deltaf(i))), zeros(1,4)];
    alphaFRLowerBound = [
        (alphaMin*sin(deltaf(i)) - cos(deltaf(i))), (sin(deltaf(i)) + alphaMin*cos(deltaf(i))), 0,...
        (alphaMin*(w*cos(deltaf(i)) + lf*sin(deltaf(i))) - lf*cos(deltaf(i)) + w*sin(deltaf(i))), zeros(1,4)];
    stateConstraintBlock = [
        slipFLUpperBound; slipFRUpperBound; slipRLUpperBound; slipRRUpperBound;...
        slipFLLowerBound; slipFRLowerBound; slipRLLowerBound; slipRRLowerBound;...
        alphaFLUpperBound; alphaFRUpperBound; alphaRLUpperBound; alphaRRUpperBound;...
        alphaFLLowerBound; alphaFRLowerBound; alphaRLLowerBound; alphaRRLowerBound];
    stateConstraintBlock = stateConstraintBlock * Mxi;
    AStateConstraint = [kron(speye(N+1), stateConstraintBlock), zeros(neps*(N+1),2*nu*M), -speye(neps*(N+1))];
    
    % update the problem matrices/vectors
    A = []; l = []; u = [];
    if constraints(1); A = [A; ADynamics]; l = [l; lDynamics]; u = [u; uDynamics]; end
    if constraints(2); A = [A; AStateConstraint]; l = [l; lStateConstraint]; u = [u; uStateConstraint]; end  
    if constraints(3); A = [A; AInputConstraint]; l = [l; lInputConstraint]; u = [u; uInputConstraint]; end
    if constraints(4); A = [A; AInputRateConstraint]; l = [l; lInputRateConstraint]; u = [u; uInputRateConstraint]; end
    if constraints(4); A = [A; AdeltaU]; l = [l; ldeltaU]; u = [u; udeltaU]; end
    if constraints(2); A = [A; ASlack]; l = [l; lSlack]; u = [u; uSlack]; end
    
    % OSQP expects a vector of values for updating A; we need zeros too
    prob.update('q', q, 'Ax', full(A(Ax_idx)), 'l', l, 'u', u);
    
    % Solve
    res = prob.solve();

    % Check solver status
    if ~strcmp(res.info.status, 'solved')
        error(['OSQP did not solve the problem! ->  ',res.info.status])
    end
    
    predictions = reshape(res.x(1:nx*(N+1)),nx,N+1);
%     % plot predictions
%     cla;
%     stairs(i:i+N,xiref(2,i:i+N),'k')
%     stairs(i:i+N,predictions(2,:),'b')
%     stairs(1:i,Xis(2,1:i),'g');title('Prediction+reference+state')
%     xlabel('k')
%     ylabel('vx [m/s')
%     legend('reference','prediction','trajectory')
%     pause(0.01)
    
    % apply the inputs to the plant
    U = res.x(nx*(N+1)+1:nx*(N+1)+nu);  % take the first optimal input
    U = Mu * U  % convert back from Pi-space
    Xi = Mxi * Xi;
    [~,Xi] = ode45(@OdeSimulationModel,[0,Ts],[Xi;U;deltaf(i);ay;ax]);
    Xi = Xi(end,1:nx)'
    i
    
    % calculate the acceleration
    f = OdeSimulationModel(0,[Xi;U;deltaf(i);ay;ax]);
    ay = f(1); ax = f(2);
    
    % save optimization results (normal space) for the next iteration
    oldXis = [Xi, Mxi * predictions(:,2:end)];
    oldUs = Mu * reshape(res.x(nx*(N+1)+1:nx*(N+1)+nu*M),nu,M);
    
    % save data for plotting
    Xis(:,i+1) = Xi;
    Us(:,i) = U;
end
toc

%% plot the results
t = 0:Ts:nsim*Ts;
xiref = Mxi * xiref(:,1:nsim+1);

% states and reference
figure('Name', 'states')
ylabels = {'vy [m/s]','vx [m/s]','$\theta$ [rad]','$\dot{\theta}$ [rad/s]',...
           '$\omega_{FL}$ [rad/s]','$\omega_{FR}$ [rad/s]','$\omega_{RL}$ [rad/s]','$\omega_{RR}$ [rad/s]'};
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

% % wheel speeds
% figure
% for i=5:8
%     subplot(4,1,i-4)
%     plot(t,Xis(i,:))
%     xlabel('Time [s]')
%     ylabel(ylabels{i})
% end

% slip and slip angle
figure('Name', 'slips')
subplot(2,1,1)
sx = params.rw*Xis(5:8,:)./Xis(2,:)-1;
plot(t,100*sx')
xlabel('Time [s]')
ylabel('Longitudinal slip [\%]')
title('Longitudinal slip')
legend('\kappa_{fl}','\kappa_{fr}','\kappa_{rl}','\kappa_{rr}')
subplot(2,1,2)
slipAngleFront = (Xis(1,:)+params.lf*Xis(4,:))./Xis(2,:)-[0 deltaf];
slipAngleRear = (Xis(1,:)-params.lf*Xis(4,:))./Xis(2,:);
plot(t,rad2deg([slipAngleFront; slipAngleRear])')
xlabel('Time [s]')
ylabel('Slip angle [$^\circ$]')
title('Slip angle')
legend('\alpha_{f}','\alpha_{r}')

% inputs
figure('Name', 'input')
stairs(t,[u0 Us]')
legend('T_{fl}','T_{fr}','T_{rl}','T_{rr}')
title('Input torque')
xlabel('Time [s]')
ylabel('Input torque [Nm]')

% input rate
figure('Name', 'input rate')
% if enablePi
%     dUs = 1/Tsp*diff([u0 Us],1,2);
% else
    dUs = 1/Ts*diff([u0 Us],1,2);
% end
stairs(t(2:end),dUs')
legend('\Delta T_{fl}','\Delta T_{fr}','\Delta T_{rl}','\Delta T_{rr}')
title('Input torque rate')
xlabel('Time [s]')
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

figure('Name', 'path')
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

% block = speye(nu*M) + kron(diag(ones(M-1,1),1),-speye(nu));
% block = block(1:nu*(M-1),:);
% AdeltaU = [zeros(nu,nx*(N+1)), -speye(nu), zeros(nu,nu*(M-1)), speye(nu), zeros(nu,nu*(M-1)+neps*(N+1));
%            zeros(nu*(M-1),nx*(N+1)), block, zeros(nu*(M-1),nu), speye(nu*(M-1)), zeros(nu*(M-1),neps*(N+1))];
