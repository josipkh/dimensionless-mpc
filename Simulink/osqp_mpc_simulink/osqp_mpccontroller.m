function ctrl = osqp_mpccontroller(x_current)
persistent lb ub
persistent Ad Bd

if isempty(Ad)
    % Discrete time model of a quadcopter
    Ad = [1       0       0   0   0   0   0.1     0       0    0       0       0;
        0       1       0   0   0   0   0       0.1     0    0       0       0;
        0       0       1   0   0   0   0       0       0.1  0       0       0;
        0.0488  0       0   1   0   0   0.0016  0       0    0.0992  0       0;
        0      -0.0488  0   0   1   0   0      -0.0016  0    0       0.0992  0;
        0       0       0   0   0   1   0       0       0    0       0       0.0992;
        0       0       0   0   0   0   1       0       0    0       0       0;
        0       0       0   0   0   0   0       1       0    0       0       0;
        0       0       0   0   0   0   0       0       1    0       0       0;
        0.9734  0       0   0   0   0   0.0488  0       0    0.9846  0       0;
        0      -0.9734  0   0   0   0   0      -0.0488  0    0       0.9846  0;
        0       0       0   0   0   0   0       0       0    0       0       0.9846];
end

if isempty(Bd)
    Bd = [0      -0.0726  0       0.0726;
        -0.0726  0       0.0726  0;
        -0.0152  0.0152 -0.0152  0.0152;
        0      -0.0006 -0.0000  0.0006;
        0.0006  0      -0.0006  0;
        0.0106  0.0106  0.0106  0.0106;
        0      -1.4512  0       1.4512;
        -1.4512  0       1.4512  0;
        -0.3049  0.3049 -0.3049  0.3049;
        0      -0.0236  0       0.0236;
        0.0236  0      -0.0236  0;
        0.2107  0.2107  0.2107  0.2107];
end

persistent N
if isempty(N)
    % Prediction horizon
    N = 10;
end

persistent nx nu
[nx, nu] = size(Bd);

persistent prob
if isempty(prob) 
    % Constraints
    u0 = 10.5916;
    umin = [9.6; 9.6; 9.6; 9.6] - u0;
    umax = [13; 13; 13; 13] - u0;
    xmin = [-pi/6; -pi/6; -Inf; -Inf; -Inf; -1; -Inf(6,1)];
    xmax = [ pi/6;  pi/6;  Inf;  Inf;  Inf; Inf; Inf(6,1)];
    
    % Objective function
    Q = diag([0 0 10 10 10 10 0 0 0 5 5 5]);
    QN = Q;
    R = 0.1*eye(4);
    
    % Initial and reference states
    x0 = zeros(12,1);
    xr = [0; 0; 1; 0; 0; 0; 0; 0; 0; 0; 0; 0];
    
    % Cast MPC problem to a QP: x = (x(0),x(1),...,x(N),u(0),...,u(N-1))
    % - quadratic objective
    P = blkdiag( kron(speye(N), Q), QN, kron(speye(N), R) );
    % - linear objective
    q = [repmat(-Q*xr, N, 1); -QN*xr; zeros(N*nu, 1)];
    % - linear dynamics
    Ax = kron(speye(N+1), -speye(nx)) + kron(sparse(diag(ones(N, 1), -1)), Ad);
    Bu = kron([sparse(1, N); speye(N)], Bd);
    Aeq = [Ax, Bu];
    leq = [-x0; zeros(N*nx, 1)];
    ueq = leq;
    % - input and state constraints
    Aineq = speye((N+1)*nx + N*nu);
    lineq = [repmat(xmin, N+1, 1); repmat(umin, N, 1)];
    uineq = [repmat(xmax, N+1, 1); repmat(umax, N, 1)];
    % - OSQP constraints
    A = [Aeq; Aineq];
    lb = [leq; lineq];
    ub = [ueq; uineq];
    
    % Create an OSQP object
    prob = osqp;
    
    % Setup workspace
    prob.setup(P, q, A, lb, ub, 'warm_start', true);
   
end

% Update initial state
lb(1:nx) = -x_current;
ub(1:nx) = -x_current;

prob.update('l', lb, 'u', ub);
res = prob.solve();

if ~strcmp(res.info.status, 'solved')
    error('OSQP did not solve the problem!')
end

ctrl = res.x((N+1)*nx+1:(N+1)*nx+nu);

end
