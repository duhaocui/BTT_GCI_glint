clear
clc
close all

DEMO_FLAG = 0;

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
x_R = 0;
y_R = 0;
sigma_r = 100;
sigma_theta = 0.05;
meas_param{1} = sigma_r;
meas_param{2} = sigma_theta;
meas_param{3} = x_R;
meas_param{4} = y_R;

%% UKF
beta_UKF = beta_tgt;
n_x = size(F, 1);
alpha = 1;
beta = 0;
% kappa = 3 - n_x;
kappa = 0;
tr_param{1} = alpha;
tr_param{2} = beta;
tr_param{3} = kappa;

%% main function
[mse_px, mse_py] = func_BTT_UKF(dyn_param, meas_param, tr_param, DEMO_FLAG);