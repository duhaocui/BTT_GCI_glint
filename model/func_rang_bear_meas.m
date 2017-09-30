function  z = func_rang_bear_meas(pos, meas_param)

sigma_r = meas_param{1};
sigma_theta = meas_param{2};
x_R = meas_param{3};
y_R = meas_param{4};

x = pos(1);
y = pos(2);

v_r = normrnd(0, sigma_r);
v_theta = normrnd(0, sigma_theta);
z_r = sqrt( (x - x_R)^2 + (y - y_R)^2) + v_r;
z_theta = atan2( (y - y_R), (x - x_R) ) + v_theta;

z = zeros(2, 1);
z(1) = z_r;
z(2) = z_theta;