function [x_pred, P_pred] = func_EKF_pred(x_upd, P_upd, F, G, g, f_J, Q, beta_EKF)

f = func_calc_air_dens(x_upd, g, beta_EKF);
x_pred = F * x_upd + G * f + G * [0; -g];
P_pred = (F + G * f_J) * P_upd * (F + G * f_J)' + Q;

