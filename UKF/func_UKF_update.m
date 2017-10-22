function [m_upd, P_upd] = func_UKF_update(m_pred, P_pred, z, func_handle, meas_param, tr_param)

sigma_r = meas_param{1};
sigma_theta = meas_param{2};
R = diag([sigma_r^2, sigma_theta^2]);

[m_tfr, P_tfr, P_cov_tfr] = func_ut_transform(m_pred, P_pred, func_handle, meas_param, tr_param);

% testing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [m_tfr_1, P_tfr_1, P_cov_tfr_1] = ut_transform(m_pred, P_pred, func_handle, meas_param, tr_param);
% norm(m_tfr - m_tfr_1)
% norm(P_tfr - P_tfr_1)
% norm(P_cov_tfr - P_cov_tfr_1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P_tfr = P_tfr + R;
K = P_cov_tfr / P_tfr;
m_upd = m_pred + K * (z - m_tfr);
P_upd = P_pred - K * P_tfr * K';