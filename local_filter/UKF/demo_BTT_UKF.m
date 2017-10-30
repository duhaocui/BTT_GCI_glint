clear
clc
close all

DEMO_FLAG = 1;

%% target parameters
T = 2;
q = 1;
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
dyn_param{1} = g;
dyn_param{2} = beta_tgt;
dyn_param{3} = F;
dyn_param{4} = G;
dyn_param{5} = T;
dyn_param{6} = q;

%% measurement parameters
x_s = 0;
y_s = 0;
r_range = 10;
r_theta = 0.05;
R = diag([r_range^2, r_theta^2]);
meas_param{1} = x_s;
meas_param{2} = y_s;
meas_param{3} = R;

%% unscented transform parameters
alpha = 1;
beta = 0;
% kappa = 3 - n_x;
kappa = 0;
ut_param{1} = alpha;
ut_param{2} = beta;
ut_param{3} = kappa;

%% main function
[mse_px, mse_py] = func_BTT_UKF(dyn_param, meas_param, ut_param, DEMO_FLAG);