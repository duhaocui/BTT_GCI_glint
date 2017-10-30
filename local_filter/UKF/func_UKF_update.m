function [x_upd, P_upd] = func_UKF_update(x_pred, P_pred, z, R, func_handle, meas_param, ut_param)

% construct sigma points
sigma_pred = func_constr_sigma(x_pred, P_pred, ut_param);

% read data from `sigma_pred`
X_in = sigma_pred.X;
W_m = sigma_pred.W_m;
W_c = sigma_pred.W_c;

% unscented transform
X_ut = func_ut_transform(X_in, func_handle, meas_param);

%% calculate `x_ut`, `P_ut` and `P_cov_ut`
L = size(X_in, 2);
% calculate `m_ut`
x_ut = [];
for i = 1 : L
    if i == 1
        x_ut = W_m(i) * X_ut(:, i);
    else
        x_ut = x_ut + W_m(i) * X_ut(:, i);
    end
end
% calculate `P_ut`
P_ut = [];
for i = 1 : L
    if i == 1
        P_ut = W_c(i) * (X_ut(:, i) - x_ut) * (X_ut(:, i) - x_ut)';
    else
        P_ut = P_ut + W_c(i) * (X_ut(:, i) - x_ut) * (X_ut(:, i) - x_ut)';
    end
end
% calculate `P_cov_ut`
P_cov_ut = [];
for i = 1 : L
    if i == 1
        P_cov_ut = W_c(i) * (X_in(:, i) - x_pred) * (X_ut(:, i) - x_ut)';
    else
        P_cov_ut = P_cov_ut + W_c(i) * (X_in(:, i) - x_pred) * (X_ut(:, i) - x_ut)';
    end
end

%% calcualte `m_upd` and `P_upd`
P_ut = P_ut + R;
K = P_cov_ut / P_ut;
x_upd = x_pred + K * (z - x_ut);
P_upd = P_pred - K * P_ut * K';