function [m_pred, P_pred] = func_UKF_predict(m, P, func_handle, dyn_param, tr_param)

T = dyn_param{5};
q = dyn_param{6};
Theta = [T^3 / 3, T^2 / 2; T^2 / 2, T];
Q = q * blkdiag(Theta, Theta);

[m_tfr, P_tfr] = func_ut_transform(m, P, func_handle, dyn_param, tr_param);

% testing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [m_tfr_1, P_tfr_1] = ut_transform(m, P, func_handle, dyn_param, tr_param);
% norm(m_tfr - m_tfr_1)
% norm(P_tfr - P_tfr_1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m_pred = m_tfr;
P_pred = P_tfr + Q;

