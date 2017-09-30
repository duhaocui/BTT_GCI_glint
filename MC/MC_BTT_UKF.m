clear
clc
close all

cd('/Users/Colin/Documents/GitHub/BTT_GCI_glint')

DEMO_FLAG = 0;
mcruns = 1000;

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
n = size(F, 1);
alpha = 1;
beta = 0;
kappa = 3 - n;
tr_param{1} = alpha;
tr_param{2} = beta;
tr_param{3} = kappa;

%% Monte-carlo simulation
for i = 1 : mcruns
    disp(sprintf('Monte Carlo runs: NO. %d', i) );
    [mse_px, mse_py] = func_BTT_UKF(T, dyn_param, q, meas_param, sigma_r, sigma_theta, tr_param, beta_UKF, DEMO_FLAG);
    if i == 1
        MSE_px = mse_px;
        MSE_py = mse_py;
    else
        MSE_px = MSE_px + mse_px;
        MSE_py = MSE_py + mse_py;
    end
end

MSE_px = MSE_px / mcruns;
MSE_py = MSE_py / mcruns;

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

% log
logTime = sprintf('%s', datestr(now,30));
if ismac
    dataName = strcat('log/data_UKF_', logTime);
else
    dataName = strcat('log\data_UKF_', logTime);
end
save(dataName);