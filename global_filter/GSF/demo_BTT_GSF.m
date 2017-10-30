clear
clc
close all

DEMO_FLAG = 0;

%% target parameters
T = 2;
q = [1; 0.01];
g = 9.81;
beta_tgt = 4 * 1e4;
F = [1, T, 0, 0;
    0, 1, 0, 0;
    0, 0, 1, T;
    0, 0, 0, 1];
G = [T^2 / 2, 0
    T,       0;
    0,       T^2 / 2;
    0,       T];
epsilon = 0.5;
dyn_param{1} = g;
dyn_param{2} = beta_tgt;
dyn_param{3} = F;
dyn_param{4} = G;
dyn_param{5} = T;
dyn_param{6} = q;
assert(length(q) == 2)
dyn_param{7} = epsilon;

%% measurement parameters
x_R = 0;
y_R = 0;
sigma_r = [100; 1];
sigma_theta = [0.05; 0.005];
meas_param{1} = sigma_r;
meas_param{2} = sigma_theta;
meas_param{3} = x_R;
meas_param{4} = y_R;
assert(length(sigma_r) == 2)
assert(length(sigma_theta) == 2)
meas_param{5} = epsilon;