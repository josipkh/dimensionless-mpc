%% Bicycle model & Pi groups - MPC controller
clear;clc;close all;

%% load vehicle parameters
vehicle_params = 2;  % 1: lab scale, 2: full size
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

%% MPC setup
yalmip('clear');

% MPC setup
pigroups = 1;  % binary

nx = size(A,1);  % no. of states
nu = size(B,2);  % no. of inputs
Cref = [1 0 0 0];  % tracked states (only lateral position)
% if pigroups; Crefp = Cref*M; end
ny = size(Cref,1);  % no. of tracked states

Q = 10^6*diag([1 0 0 0]);  % state weights
if pigroups; Q = M*Q*M; end
QN = Q;  % terminal cost
R = 1e-6;  % input weights
S = 1e-6;  % input rate weights
N = 50;  % prediction horizon
Nu = 5;  % control horizon

% YALMIP data
x = sdpvar(nx*ones(1,N+1),ones(1,N+1));     % state variables
r = sdpvar(ny*ones(1,N+1),ones(1,N+1));     % reference variables
u = sdpvar(nu*ones(1,Nu),ones(1,Nu));       % input variables
pastu = sdpvar(1);                          % previous input
deltau = diff([pastu u{:}]);                % input rate
ops = sdpsettings('verbose',0,'solver','osqp');  % print output

% constraints and cost
constraints = [];
% steering rate limit
steering_rate_limit = VEHICLE.MAX_STEERING_RATE * Ts / VEHICLE.STEERING_RATIO;
constraints = [-steering_rate_limit <= deltau <= steering_rate_limit];
% steering angle limit
steering_angle_limit = VEHICLE.MAX_STEERING_ANGLE / VEHICLE.STEERING_RATIO;
constraints = [constraints, -steering_angle_limit <= [u{:}]<= steering_angle_limit];
objective = 0;
for k = 1:N+1
    if k == N+1
        % terminal cost
        objective = objective + (Cref*x{k}-r{k})'*QN(1)*(Cref*x{k}-r{k});
    else
        if k > Nu
            % input cost
            objective = objective + u{Nu}'*R*u{Nu};
            % state constraints
            if pigroups
                constraints = [constraints, x{k+1} == Apd*x{k} + Bpd*u{Nu}];
            else
                constraints = [constraints, x{k+1} == Ad*x{k} + Bd*u{Nu}];
            end
        else
            % input cost
            objective = objective + u{k}*R*u{k};
            % input rate cost
            objective = objective + deltau(k)*S*deltau(k);
            % state constraints
            if pigroups
                constraints = [constraints, x{k+1} == Apd*x{k} + Bpd*u{k}];
            else
                constraints = [constraints, x{k+1} == Ad*x{k} + Bd*u{k}];
            end
        end
        % tracking cost
        objective = objective + (Cref*x{k}-r{k})'*Q(1)*(Cref*x{k}-r{k});
    end
end

solutions_out = {[u{:}], [x{:}]};
parameters_in = {x{1},[r{:}],pastu};
controller = optimizer(constraints, objective, ops, parameters_in,solutions_out);

%% simulation
clf; hold on
n = 1000;
x = [0;0;deg2rad(0);0];  % initial state
xhist = nan(nx, n+1);
xhist(:,1) = x;
uhist = nan(nu, n);
uprev = 0;

rhist = nan(ny, n);
tsim = 0:Ts:(n+N)*Ts;
wavelength = 32 * (Lf+Lr);
frequency = V / wavelength;
yref = 1*(Lf+Lr)*sin(2*pi*frequency*(tsim-Ts*n/5));
yref(1:n/5) = 0;
% plot(tsim,yref)
% df_ref = pi/2*sin(4*tsim)/15;  % SWA max +/- 90 deg, steering ratio = 15
% yaw_rate_ref = V/(Lf+Lr)*tan(df_ref);
for i = 1:n
    if pigroups; x = M\x; end  % convert init. state to pi-space
    future_r = yref(i:i+N);
%         % formulate a sinusoidal y-reference
% %         wavelength = 4 * (Lf+Lr);
% %         frequency = V / wavelength;
% %         future_r = 0.5*(Lf+Lr)*sin(2*pi*frequency/Ts*(i:i+N));
%         future_r = 0.5*(Lf+Lr)*sin((i:i+N)/20);
%         future_r = yaw_rate_ref(i)*ones(1,N+1);  % assume constant ref. yaw rate

    % convert the reference to pi-space
    if pigroups; future_r = future_r/M(1); end
    inputs = {x,future_r,uprev};
    
    % find the optimal input
    [solutions,diagnostics] = controller{inputs};    
    U = solutions{1};
    X = solutions{2};
    if diagnostics ~= 0
        error('Solver error');
    end    
    
    % live plot
%     subplot(2,1,1);stairs(i:i+length(U)-1,180/pi*U,'r');title('Predicted input [deg]')
%     subplot(2,1,2);cla;stairs(i:i+N,X(1,:),'b');
%     hold on;
%     stairs(i:i+N,future_r,'k');
%     if pigroups
%         stairs(1:i,1/M(1)*xhist(1,1:i),'g');
%     else
%         stairs(1:i,xhist(1,1:i),'g');
%     end
%     title('Prediction, reference, state')
    
    % convert back to "normal" space
    if pigroups; x = M*x; end  
    
    % simulate one time step
    if vehicle_params == 1
        [~,x]=ode45(@OdeBicycleModelSmall,[0,Ts],[x;U(1)]);
    elseif vehicle_params == 2
        [~,x]=ode45(@OdeBicycleModelLarge,[0,Ts],[x;U(1)]);
    else
        error('no simulation model')
    end
    x = x(end,1:4)';
    
    % save for plotting
    xhist(:,i+1) = x;
    uhist(:,i) = U(1);
    uprev = U(1);
    
    if pigroups; future_r = future_r*M(1); end  % convert back
    rhist(:,i) = future_r(:,1);
    
    % display progress
    if ~mod(i,50) || i==1; disp(i); end
    pause(0.01)   
end

%% plotting
figure; t = (1:n)*Ts;
subplot(5,1,1); plot(t,180/pi*uhist*VEHICLE.STEERING_RATIO); title('Steering angle [deg]')
    hold on; yline(180/pi*VEHICLE.MAX_STEERING_ANGLE,'--r'); yline(-180/pi*VEHICLE.MAX_STEERING_ANGLE,'--r')
subplot(5,1,2); plot(t,[0 180/pi*diff(uhist)/Ts*VEHICLE.STEERING_RATIO]); title('Steering angle rate [deg/s]')
    hold on; yline(180/pi*VEHICLE.MAX_STEERING_RATE,'--r'); yline(-180/pi*VEHICLE.MAX_STEERING_RATE,'--r')
subplot(5,1,3); plot([0 t],xhist(1,:)); title('Lateral position [m]')
    hold on; plot(t,rhist(1,:),'r--')
subplot(5,1,4); plot([0 t],180/pi*xhist(3,:)); title('Yaw angle [deg]')
subplot(5,1,5); plot([0 t],180/pi*xhist(4,:),'b'); title('Yaw rate [deg/s]')
xlabel('t [s]')

if pigroups
    figure
    subplot(5,1,1); plot(t,180/pi*uhist*VEHICLE.STEERING_RATIO); title('Steering angle [deg]')
        hold on; yline(180/pi*VEHICLE.MAX_STEERING_ANGLE,'--r'); yline(-180/pi*VEHICLE.MAX_STEERING_ANGLE,'--r')
    subplot(5,1,2); plot(t,[0 180/pi*diff(uhist)/Ts*VEHICLE.STEERING_RATIO]); title('Steering angle rate [deg/s]')
        hold on; yline(180/pi*VEHICLE.MAX_STEERING_RATE,'--r'); yline(-180/pi*VEHICLE.MAX_STEERING_RATE,'--r')
    subplot(5,1,3); plot([0 t],1/M(1)*xhist(1,:)); title('Lateral position [-]')
        hold on; plot(t,1/M(1)*rhist(1,:),'r--')
    subplot(5,1,4); plot([0 t],180/pi*xhist(3,:)); title('Yaw angle [deg]')
    subplot(5,1,5); plot([0 t],1/M(end)*xhist(4,:),'b'); title('Yaw rate [-]')
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
