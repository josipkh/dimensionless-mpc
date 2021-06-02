%% Nonlinear MPC
clear;clc
Ts = 0.05;
nlobj = nlmpc(8,8,4);
nlobj.Ts = Ts;
nlobj.PredictionHorizon = 5;
nlobj.ControlHorizon = 2;

% system dynamics
nlobj.Model.StateFcn = "twoTrackModel";
nlobj.Model.IsContinuousTime = true;
nlobj.Jacobian.StateFcn = "twoTrackModelJacobian";
%nlobj.Model.OutputFcn = @(x,u) [x(1);x(2);x(3);x(4)];
%nlobj.Jacobian.OutputFcn = @(x,u) [eye(4) zeros(4)];

% cost function weights
nlobj.Weights.OutputVariables = [10 100 1e4 2e4 0 0 0 0];
nlobj.Weights.ManipulatedVariables = 1e-3*ones(1,4);
nlobj.Weights.ManipulatedVariablesRate = 1e-2*ones(1,4);

% constraints
for ct = 1:4
    nlobj.MV(ct).Min = -500;
    nlobj.MV(ct).Max = 500;
    nlobj.MV(ct).RateMin = -50;
    nlobj.MV(ct).RateMax = 50;
end
% limit tire slip and slip angle
% create a custom constraint function?
% nlobj.OV(1).Min = -10;
% nlobj.OV(1).Max = 10;

% set scale factors
% nlobj.OV(1).ScaleFactor = 15;   % Typical value of longitudinal velocity
% nlobj.OV(2).ScaleFactor = 0.5;  % Range for lateral deviation
% nlobj.OV(3).ScaleFactor = 0.5;  % Range for relative yaw angle
% nlobj.MV(1).ScaleFactor = 6;    % Range of steering angle
% nlobj.MV(2).ScaleFactor = 2.26; % Range of acceleration
% nlobj.MD(1).ScaleFactor = 0.2;  % Range of Curvature

x0 = [5 5 0 0 10 10 10 10]';
u0 = [0 0 0 0]';
validateFcns(nlobj,x0,u0);

%%
clc
x = x0;
u = u0;
Duration = 5;
xHistory = x;
for ct = 1:(Duration/Ts)
    % Compute optimal control moves
    [u,nloptions] = nlmpcmove(nlobj,x,u);
    % Implement first optimal control move
    x = x+Ts*twoTrackModel(x,u);    
    % Save plant states
    xHistory = [xHistory x];
    ct*Ts
end