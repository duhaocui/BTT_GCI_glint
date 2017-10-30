function [x_upd, P_upd] = func_EKF_upd(x_R, y_R, x_pred, P_pred, z, h_J, R)

x = x_pred(1);
y = x_pred(3);

I = eye(length(x_pred) );
K = P_pred * h_J' * inv(h_J * P_pred * h_J' + R);

z_r = sqrt( (x - x_R)^2 + (y - y_R)^2);
z_theta = atan2( (y - y_R), (x - x_R) );

x_upd = x_pred + K * (z - [z_r; z_theta]);
P_upd = (I - K * h_J) * P_pred;