clear
clc
close all

cd('/Users/admin/Documents/GitHub/BTT_GCI_glint')

DEMO_FLAG = 0;
mcruns = 2;

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

%% Monte-Carlo simulation
for i = 1 : mcruns
    disp(sprintf('Monte Carlo runs: NO. %d', i) );
    [mse1_px, mse1_py, beta1_est, mse2_px, mse2_py, beta2_est, mse_f_px, mse_f_py] = func_BTT_GM_GCI_EA(T, q, g, beta_tgt, F, G, nSensors, p_R, sigma_r, sigma_theta, nModes, beta1_IMM, beta2_IMM, DEMO_FLAG);
    if i == 1
        MSE1_px = mse1_px;
        MSE2_px = mse2_px;
        MSE_f_px = mse_f_px;
        MSE1_py = mse1_py;
        MSE2_py = mse2_py;
        MSE_f_py = mse_f_py;
        BETA1_est = beta1_est;
        BETA2_est = beta2_est;
    else
        MSE1_px = MSE1_px + mse1_px;
        MSE2_px = MSE2_px + mse2_px;
        MSE_f_px = MSE_f_px + mse_f_px;
        MSE1_py = MSE1_py + mse1_py;
        MSE2_py = MSE2_py + mse2_py;
        MSE_f_py = MSE_f_py + mse_f_py;
        BETA1_est = BETA1_est + beta1_est;
        BETA2_est = BETA2_est + beta2_est;
    end
end

MSE1_px = MSE1_px / mcruns;
MSE2_px = MSE2_px / mcruns;
MSE_f_px = MSE_f_px / mcruns;
MSE1_py = MSE1_py / mcruns;
MSE2_py = MSE2_py / mcruns;
MSE_f_py = MSE_f_py / mcruns;
BETA1_est = BETA1_est / mcruns;
BETA2_est = BETA2_est / mcruns;

RMSE1_px = sqrt(MSE1_px);
RMSE2_px = sqrt(MSE2_px);
RMSE_f_px = sqrt(MSE_f_px);
RMSE1_py = sqrt(MSE1_py);
RMSE2_py = sqrt(MSE2_py);
RMSE_f_py = sqrt(MSE_f_py);

figure
subplot(311)
plot(T * (1 : length(RMSE1_px) ), RMSE1_px)
xlabel('time (s)')
ylabel('RMSE1 in X (m)')
grid on
subplot(312)
plot(T * (1 : length(RMSE2_px) ), RMSE2_px)
xlabel('time (s)')
ylabel('RMSE2 in X (m)')
grid on
subplot(313)
plot(T * (1 : length(RMSE_f_px) ), RMSE_f_px)
xlabel('time (s)')
ylabel('RMSE_f in X (m)')
grid on

figure
subplot(311)
plot(T * (1 : length(RMSE1_py) ), RMSE1_py)
xlabel('time (s)')
grid on
ylabel('RMSE in Y (m)')
subplot(312)
plot(T * (1 : length(RMSE2_py) ), RMSE2_py)
xlabel('time (s)')
grid on
ylabel('RMSE in Y (m)')
subplot(313)
plot(T * (1 : length(RMSE_f_py) ), RMSE_f_py)
xlabel('time (s)')
ylabel('RMSE_f in Y (m)')
grid on

figure
subplot(211)
plot(T * (1 : length(BETA1_est) ), BETA1_est)
xlabel('time (s)')
ylabel('average values of the estimated ballistic coefficient')
grid on
subplot(212)
plot(T * (1 : length(BETA2_est) ), BETA2_est)
xlabel('time (s)')
ylabel('average values of the estimated ballistic coefficient')
grid on