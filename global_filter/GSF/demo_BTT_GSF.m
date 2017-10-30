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
epsilon_w = 0.5;
dyn_param{1} = g;
dyn_param{2} = beta_tgt;
dyn_param{3} = F;
dyn_param{4} = G;
dyn_param{5} = T;
dyn_param{6} = q;
assert(length(q) == 2)
dyn_param{7} = epsilon_w;

%% measurement parameters
% sensor position (x_s, y_s)
x_s = 0;
y_s = 0;
r_range = [10; 1];
r_theta = [0.05; 0.005];
assert(length(r_range) == 2)
assert(length(r_theta) == 2)
R = cell(2, 1);
R{1} = diag([(r_range(1) )^2, (r_theta(1) )^2]);
R{2} = diag([(r_range(2) )^2, (r_theta(2) )^2]);
epsilon_v = 0.5;
meas_param{1} = x_s;
meas_param{2} = y_s;
meas_param{3} = R;
meas_param{4} = epsilon_v;

%% unscented transform parameters
alpha = 1;
beta = 0;
% kappa = 3 - n_x;
kappa = 0;
ut_param{1} = alpha;
ut_param{2} = beta;
ut_param{3} = kappa;

%% main function
[mse_px, mse_py] = func_BTT_GSF(dyn_param, meas_param, ut_param, DEMO_FLAG);