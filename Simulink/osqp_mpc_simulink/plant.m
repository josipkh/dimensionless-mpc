clear
close all
clc

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