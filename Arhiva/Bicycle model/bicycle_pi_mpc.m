%% Bicycle model (2001) & Pi groups - MPC controller
clear;clc;

%% load vehicle parameters
vehicle_params = 1;  % 1: lab-scale IRS(2001), 2: full size (Carmaker, demo compact)
f = 1;  % scaling factor
load_vehicle_parameters;

%% system dynamics
if vehicle_params == 1
    Ts = 0.01;  % sampling time
elseif vehicle_params == 2
    Ts = 0.01;  % sampling time
end
formulate_dynamics;  % A,B,C,Ad,Bd

% Pi-groups
formulate_pi_groups;  % M,Ap,Bp,Cp,Apd,Bpd

%% MPC
yalmip('clear');close all;clc

% MPC setup
reference = 1;  % binary
pigroups = 1;  % binary

nx = size(A,1);  % no. of states
nu = size(B,2);  % no. of inputs
Cref = [0 0 0 1];  % tracked states (only yaw rate)
if pigroups; Crefp = Cref*M; end
ny = size(Cref,1);  % no. of tracked states

Q = 10^6*diag([0 0 0 1]);  % state weights
if pigroups; Q = M*Q*M; end
R = 1e-6;  % input weights
N = 20;  % horizon length

% YALMIP data
u = sdpvar(nu*ones(1,N),ones(1,N));
x = sdpvar(nx*ones(1,N+1),ones(1,N+1));
r = sdpvar(ny*ones(1,N+1),ones(1,N+1));
ops = sdpsettings('verbose',0,'solver','osqp');  % print output

% constraints and cost
constraints = [];
pastu = sdpvar(1);
constraints = [-6.28/15 <= diff([pastu u{:}]) <= 6.28/15];  % constrain du, steering rate limit
objective = 0;
for k = 1:N
    if reference
        objective = objective + (Cref*x{k}-r{k})'*Q(end)*(Cref*x{k}-r{k}) + u{k}'*R*u{k};
    else
        objective = objective + x{k}'*Q*x{k} + u{k}'*R*u{k};
    end
    if pigroups
        constraints = [constraints, x{k+1} == Apd*x{k} + Bpd*u{k}];
    else
        constraints = [constraints, x{k+1} == Ad*x{k} + Bd*u{k}];
    end
    constraints = [constraints, -6.28/15 <= u{k}<= 6.28/15];  % steering angle limit
end
if reference
    objective = objective + (Cref*x{N+1}-r{N+1})'*(Cref*x{N+1}-r{N+1});
else
    objective = objective + x{N+1}'*Q*x{N+1};
end

solutions_out = {[u{:}], [x{:}]};
if reference
    parameters_in = {x{1},[r{:}]};
    x = [0;0;0;0];  % initial state
else
    parameters_in = x{1};
    if vehicle_params == 1
        x = [0.1;0;0.1;0];  % initial state, lab-scale vehicle
    else
        x = [0.2;0;0;0];  % full size vehicle
    end    
end

controller = optimizer(constraints, objective, ops, parameters_in,solutions_out);

clf; hold on
n = 250;
xhist = nan(nx, n+1);
xhist(:,1) = x;
uhist = nan(nu, n);
if reference; rhist = nan(ny, n); end
% wavelength = 4 * (Lf+Lr);
% frequency = V / wavelength;
% yref = 0.5*(Lf+Lr)*sin(2*pi*frequency*(0:0.01:n*Ts));
% figure;plot(0:0.01:n*Ts,yref)
tsim = 0:Ts:(n+N)*Ts;
df_ref = pi/2*sin(4*tsim)/15;  % SWA max +/- 90 deg, steering ratio = 15
yaw_rate_ref = V/(Lf+Lr)*tan(df_ref);
for i = 1:n
    if pigroups; x = M\x; end  % convert init. state to pi-space
    if reference
%         % formulate a sinusoidal y-reference
% %         wavelength = 4 * (Lf+Lr);
% %         frequency = V / wavelength;
% %         future_r = 0.5*(Lf+Lr)*sin(2*pi*frequency/Ts*(i:i+N));
%         future_r = 0.5*(Lf+Lr)*sin((i:i+N)/20);

        future_r = yaw_rate_ref(i)*ones(1,N+1);  % assume constant ref. yaw rate
        % convert the reference to pi-space
        if pigroups; future_r = future_r/M(end); end
        inputs = {x,future_r};
    else
        inputs = {x};
    end
    [solutions,diagnostics] = controller{inputs};    
    U = solutions{1};
    X = solutions{2};
    if diagnostics ~= 0
        error('Solver error');
    end    
    
    subplot(2,1,1);stairs(i:i+length(U)-1,180/pi*U,'r');title('Predicted input [deg]')
    subplot(2,1,2);cla;stairs(i:i+N,X(4,:),'b');
    hold on;
    if reference; stairs(i:i+N,future_r,'k'); end
    if pigroups
        stairs(1:i,1/M(end)*xhist(4,1:i),'g');
    else
        stairs(1:i,xhist(4,1:i),'g');
    end
    title('Prediction, reference, state')
    
    if pigroups; x = M*x; end  % convert back to "normal" space
    if vehicle_params == 1
        [~,x]=ode45(@OdeBicycleModelSmall,[0,Ts],[x;U(1)]);
    elseif vehicle_params == 2
        [~,x]=ode45(@OdeBicycleModelLarge,[0,Ts],[x;U(1)]);
    else
        error('no simulation model')
    end
    x = x(end,1:4)';
    
    xhist(:,i+1) = x;
    uhist(:,i) = U(1);
    if reference
        if pigroups; future_r = future_r*M(end); end  % convert back
        rhist(:,i) = future_r(:,1);
    end
    pause(0.01)   
end

figure; t = (1:n)*Ts;
subplot(4,1,1); plot(t,180/pi*uhist); title('Steering angle [deg]')
subplot(4,1,2); plot([0 t],xhist(1,:)); title('Lateral position [m]')
subplot(4,1,3); plot([0 t],180/pi*xhist(3,:)); title('Yaw angle [deg]')
subplot(4,1,4); plot([0 t],xhist(4,:),'b'); title('Yaw rate [rad/s]')
if reference
    hold on; plot(t,rhist(1,:),'r--')
end
%% LQR
% Sp = ss(Ap,Bp,Cp,0);
% Q = eye(4);
% R = 1;
% L=lqr(Sp,Q,R);
% Spcl=ss(Ap-Bp*L,Bp,Cp,0);
% Scl=ss(A-B*L/M,B,C,0);
% figure;step(Spcl)
% figure;step(Scl)

%% symbolic calculations
% syms m Iz Lf Lr Cf Cr V real
% A1 = -(Cf+Cr)/m;
% A2 = (Cr*Lr-Cf*Lf)/m;
% A3 = (Cr*Lr-Cf*Lf)/Iz;
% A4 = -(Cf*Lf^2+Cr*Lr^2)/Iz;
% 
% % system dynamics
% A = [0   1    0   0;
%      0 A1/V -A1 A2/V;
%      0   0    0   1;
%      0 A3/V -A3 A4/V];
% B = [0; Cf/m; 0; Lf*Cf/Iz];

%Ap==(Lf+Lr)/V*inv(M)*A*M  % (check difference, 10e-15)
%Bp==(Lf+Lr)/V*inv(M)*B