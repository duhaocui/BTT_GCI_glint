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
x_R = 0;
y_R = 0;
sigma_r = 100;
sigma_theta = 0.05;

%% IMM parameters
% nModes = 2;
% beta_IMM = zeros(nModes, 1);
% beta_IMM(2) = 1 * 1e4;
% beta_IMM(1) = 6 * 1e4;

nModes = 3;
beta_IMM = zeros(nModes, 1);
beta_IMM(1) = 1 * 1e4;
beta_IMM(2) = 3 * 1e4;
beta_IMM(3) = 5 * 1e4;

%% main function
[mse_px, mse_py, beta_est] = func_BTT_IMM(T, q, g, beta_tgt, F, G, x_R, y_R, sigma_r, sigma_theta, nModes, beta_IMM, DEMO_FLAG);

