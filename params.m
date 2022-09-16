clc;clear;

% Environment
Je = 0;
de = 0;
Ke = 25;

% Actuator
Jm = 0.0104;
dm = 0.0;
n = 102.8393;

% Spring
K = 3;

% Controller parameters
Kp = 20;
Ki = 0;
Kd = 0.8;

% DOB filter
wc = 100;
w2 = wc*wc;
