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

%% measurement parameters
x_R = 0;
y_R = 0;
sigma_r = 100;
sigma_theta = 0.05;

%% EKF
beta_EKF = beta_tgt;

%% main function
[mse_px, mse_py] = func_BTT_EKF(T, q, g, beta_tgt, F, G, x_R, y_R, sigma_r, sigma_theta, beta_EKF, DEMO_FLAG);