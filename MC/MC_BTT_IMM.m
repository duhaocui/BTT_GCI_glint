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
x_R = 0;
y_R = 0;
sigma_r = 100;
sigma_theta = 0.05;

%% IMM parameters
nModes = 2;
beta_IMM = zeros(nModes, 1);
beta_IMM(1) = beta_tgt;
beta_IMM(2) = 1 * 1e4;

%% Monte-carlo simulation
for i = 1 : mcruns
    disp(sprintf('Monte Carlo runs: NO. %d', i) );
    [mse_px, mse_py, beta_est] = func_BTT_IMM(T, q, g, beta_tgt, F, G, x_R, y_R, sigma_r, sigma_theta, nModes, beta_IMM, DEMO_FLAG);
    if i == 1
        MSE_px = mse_px;
        MSE_py = mse_py;
        BETA_est = beta_est;
    else
        MSE_px = MSE_px + mse_px;
        MSE_py = MSE_py + mse_py;
        BETA_est = BETA_est + beta_est;
    end
end

MSE_px = MSE_px / mcruns;
MSE_py = MSE_py / mcruns;
BETA_est = BETA_est / mcruns;

RMSE_px = sqrt(MSE_px);
RMSE_py = sqrt(MSE_py);

figure
plot(T * (1 : length(RMSE_px) ), RMSE_px)
xlabel('time (s)')
ylabel('RMSE in X (m)')
grid on
figure
plot(T * (1 : length(RMSE_py) ), RMSE_py)
xlabel('time (s)')
grid on
ylabel('RMSE in Y (m)')
figure
plot(T * (1 : length(BETA_est) ), BETA_est)
xlabel('time (s)')
ylabel('average values of the estimated ballistic coefficient')
grid on

% log
logTime = sprintf('%s', datestr(now,30));
if ismac
    dataName = strcat('log/data_IMM_', logTime);
else
    dataName = strcat('log\data_IMM_', logTime);
end
save(dataName);