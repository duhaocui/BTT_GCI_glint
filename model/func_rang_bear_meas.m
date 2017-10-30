function  z = func_rang_bear_meas(pos_tgt, meas_param)

x = pos_tgt(1);
y = pos_tgt(2);

x_s = meas_param{1};
y_s = meas_param{2};

z_range = sqrt( (x - x_s)^2 + (y - y_s)^2);
z_theta = atan2( (y - y_s), (x - x_s) );

z = zeros(2, 1);
z(1) = z_range;
z(2) = z_theta;