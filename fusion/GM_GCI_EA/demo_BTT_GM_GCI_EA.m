clear
clc
close all
DEMO_FLAG = 1;

%% target parameters
T = 1;
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

%% measurement parameters
nSensors = 2;

% [x1_R, y1_R;
%  x2_R, y2_R]
p_R = zeros(nSensors, 2);

% [sigma1_r, sigma2_r]
sigma_r = [100, 100];

% [sigma1_theta, sigma2_theta]
sigma_theta = [0.05, 0.05];

%% IMM parameters
nModes = [2, 3];

beta1_IMM = zeros(nModes(1), 1);
beta2_IMM = zeros(nModes(2), 1);

if nModes(1) == 2
    beta1_IMM(1) = 1 * 1e4;
    beta1_IMM(2) = 4 * 1e4;
elseif nModes(1) == 3
    beta1_IMM(1) = 1 * 1e4;
    beta1_IMM(2) = 3 * 1e4;
    beta1_IMM(3) = 5 * 1e4;
end

if nModes(2) == 2
    beta2_IMM(1) = 1 * 1e4;
    beta2_IMM(2) = 4 * 1e4;
elseif nModes(2) == 3
    beta2_IMM(1) = 1 * 1e4;
    beta2_IMM(2) = 3 * 1e4;
    beta2_IMM(3) = 5 * 1e4;
end

%% main function
[mse1_px, mse1_py, beta1_est, mse2_px, mse2_py, beta2_est, mse_f_px, mse_f_py] = func_BTT_GM_GCI_EA(T, q, g, beta_tgt, F, G, nSensors, p_R, sigma_r, sigma_theta, nModes, beta1_IMM, beta2_IMM, DEMO_FLAG);
