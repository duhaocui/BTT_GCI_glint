function [m_upd, P_upd] = func_UKF_update(m_pred, P_pred, z, func_handle, func_param, tr_param, R)

[m_tfr, P_tfr, P_cov_tfr] = func_ut_transform(m_pred, P_pred, func_handle, func_param, tr_param);
% [m_tfr, P_tfr, P_cov_tfr] = ut_transform(m_pred, P_pred, func_handle, func_param, tr_param);

P_tfr = P_tfr + R;
K = P_cov_tfr / P_tfr;
m_upd = m_pred + K * (z - m_tfr);
P_upd = P_pred - K * P_tfr * K';