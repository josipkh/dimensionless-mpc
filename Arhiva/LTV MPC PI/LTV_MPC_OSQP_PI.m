%% Two-track model LTV MPC
clear;clc;close all
load params.mat  % vehicle parameters
nx = 8;  % no. of states
nu = 4;  % no. of inputs
Ts = 0.05;  % sampling time
N = 10;  % horizon length
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

C = [eye(4),zeros(4)];
Cp = C * Mxi;
ny = size(C,1);  % no. of tracked states

% Constraints
u0 = zeros(nu,1);
umin = -500*ones(nu,1);
umax = 500*ones(nu,1);
% umin = -inf(nu,1);
% umax = inf (nu,1);
xmin = -inf(nx,1);
xmax = inf(nx,1);
% dumin = -50; dumax = 50

% Objective function
Q = diag([10 100 1e4 2e4 0 0 0 0]);  % tracked state weights
% Q = Mxi(1:4,1:4) * Q * Mxi(1:4,1:4);
QN = Q;  % terminal cost weights
R = 1e-3*eye(nu);  % input weights
% R = Mu * R * Mu;

% Initial and reference states
vx0 = 10;  % m/s
vx1 = 8;
xi0 = [0; vx0; 0; 0; 250/79*vx0*ones(4,1)];
xiref = repmat(xi0,1,nsim+N);  % constant forward driving
xiref = [zeros(1,nsim+N); 
         [vx0*ones(1,nsim/3),vx0+(vx1-vx0)/nsim*3*(1:nsim/3),vx1*ones(1,nsim/3+N)];
         zeros(6,nsim+N)];
% xiref = [zeros(1,nsim+N); 
%          [vx0*ones(1,20),vx1*ones(1,nsim+N-20)];
%          zeros(2,nsim+N)];
% xiref = C * Mxinv * C' * xiref;
u0 = zeros(4,1);
deltaf = zeros(nsim+N,1);

% log data
Xis = nan(nx,nsim+1);
Us = nan(nu,nsim);
Xis(:,1) = xi0;

% Simulate in closed loop
Xi = xi0;
U = u0;
figure; hold on
tic
for i = 1 : nsim
    % Linearized dynamics
    [A,B] = twoTrackModelJacobian(Xi,U,deltaf(i));
    [Ad,Bd] = c2d(A,B,Ts);
    
    % calculate predictions with constant input
    [~,xihat] = ode45(@odeTwoTrackModel,0:Ts:N*Ts,[Xi;U;deltaf(i)]);
    xihat = xihat(:,1:8)';  % nx*(N+1)
    
    % calculate linearization offset
    d = circshift(xihat,[0,-1]) - Ad*xihat - Bd*U;
    d = d(:,1:end-1);
    
    % Cast MPC problem to a QP: x = (x(0),x(1),...,x(N),u(0),...,u(N-1))
    % quadratic objective
    P = blkdiag( kron(speye(N), Q), QN, kron(speye(N), R) );
    % linear objective
    q = [reshape(-Q*xiref(:,i:i+N-1), N*nx, 1); -QN*xiref(:,i+N); zeros(N*nu, 1)];
    
    % input and state constraints
    Aineq = speye((N+1)*nx + N*nu);
    lineq = [repmat(xmin, N+1, 1); repmat(umin, N, 1)];
    uineq = [repmat(xmax, N+1, 1); repmat(umax, N, 1)];
    %Aineq = []; lineq = []; uineq = [];
    
    % system dynamics
    Ax = kron(speye(N+1), -speye(nx)) + kron(sparse(diag(ones(N, 1), -1)), Ad);
    Bu = kron([sparse(1, N); speye(N)], Bd);
    Aeq = [Ax, Bu];
    leq = [-Xi; -reshape(d, N*nx, 1)];
    ueq = leq;
    
    % OSQP constraints
    Aopt = [Aeq; Aineq];
    lopt = [leq; lineq];
    uopt = [ueq; uineq];
    
    % Create an OSQP object
    prob = osqp;

    % Setup workspace
    prob.setup(P, q, Aopt, lopt, uopt,...
               'warm_start', true,...
               'eps_abs',1e-6,...
               'eps_rel',1e-6);
    
    % Solve
    res = prob.solve();

    % Check solver status
    if ~strcmp(res.info.status, 'solved')
        error('OSQP did not solve the problem!')
    end
    
    % plot predictions
    predictions = reshape(res.x(1:(N+1)*nx),nx,N+1);
%     cla;
%     stairs(i:i+N,xiref(2,i:i+N),'k')
%     stairs(i:i+N,predictions(2,:),'b')
%     stairs(1:i,Xis(2,1:i),'g');title('Prediction+reference+state')
%     xlabel('k')
%     ylabel('vx [m/s')
%     legend('reference','prediction','trajectory')
%     pause(0.01)
    
    % apply the inputs to the plant
    U = res.x((N+1)*nx+1:(N+1)*nx+nu);
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
plot(t,[u0 Us])
legend('T_{fl}','T_{fr}','T_{rl}','T_{frr}')
xlabel('t [s]')
ylabel('Input torque [N]')