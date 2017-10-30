clear
clc
close all

DEMO_FLAG = 0;

%% target parameters
T = 2;
Theta = [T^3 / 3, T^2 / 2; T^2 / 2, T];
q = [1; 0.01];
Q_cell = cell(2, 1);
Q_cell{1} = q(1) * blkdiag(Theta, Theta);
Q_cell{2} = q(2) * blkdiag(Theta, Theta);
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
epsilon_w = 0.4;
dyn_param{1} = g;
dyn_param{2} = beta_tgt;
dyn_param{3} = F;
dyn_param{4} = G;
assert(length(q) == 2)

%% measurement parameters
% sensor position (x_s, y_s)
x_s = 0;
y_s = 0;
r_range = [10; 1];
r_theta = [0.05; 0.005];
R_cell = cell(2, 1);
R_cell{1} = diag([(r_range(1) )^2, (r_theta(1) )^2]);
R_cell{2} = diag([(r_range(2) )^2, (r_theta(2) )^2]);
epsilon_v = 0.4;
meas_param{1} = x_s;
meas_param{2} = y_s;
assert(length(r_range) == 2)
assert(length(r_theta) == 2)

%% unscented transform parameters
alpha = 1;
beta = 0;
% kappa = 3 - n_x;
kappa = 0;
ut_param{1} = alpha;
ut_param{2} = beta;
ut_param{3} = kappa;

%% main function
[mse_px, mse_py] = func_BTT_GSF(T, Q_cell, R_cell, epsilon_w, epsilon_v, dyn_param, meas_param, ut_param, DEMO_FLAG);