function  z = func_rang_bear_meas(pos, meas_param)

x_R = meas_param{1};
y_R = meas_param{2};

x = pos(1);
y = pos(2);

z_r = sqrt( (x - x_R)^2 + (y - y_R)^2);
z_theta = atan2( (y - y_R), (x - x_R) );

z = zeros(2, 1);
z(1) = z_r;
z(2) = z_theta;