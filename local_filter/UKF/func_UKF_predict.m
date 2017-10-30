function [m_pred, P_pred] = func_UKF_predict(m_upd, P_upd, Q, func_handle, dyn_param, ut_param)

% construct sigma points
sigma_upd = func_constr_sigma(m_upd, P_upd, ut_param);

% read data from `sigma_upd`
X_in = sigma_upd.X;
W_m = sigma_upd.W_m;
W_c = sigma_upd.W_c;

% unscented transform
X_ut = func_ut_transform(X_in, func_handle, dyn_param);

%% calculate `m_ut` and `P_ut`
L = size(X_in, 2);
% calculate `m_ut`
m_ut = [];
for i = 1 : L
    if i == 1
        m_ut = W_m(i) * X_ut(:, i);
    else
        m_ut = m_ut + W_m(i) * X_ut(:, i);
    end
end
% calculate `P_ut`
P_ut = [];
for i = 1 : L
    if i == 1
        P_ut = W_c(i) * (X_ut(:, i) - m_ut) * (X_ut(:, i) - m_ut)';
    else
        P_ut = P_ut + W_c(i) * (X_ut(:, i) - m_ut) * (X_ut(:, i) - m_ut)';
    end
end

%% calculate `m_pred` and `P_pred`
m_pred = m_ut;
P_pred = P_ut + Q;
