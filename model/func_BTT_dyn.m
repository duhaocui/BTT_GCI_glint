function x_curr = func_BTT_dyn(x_prev, dyn_param)

g = dyn_param{1};
beta_tgt = dyn_param{2};
q = dyn_param{3};
F = dyn_param{4};
G = dyn_param{5};
T = dyn_param{6};

Theta = [T^3 / 3, T^2 / 2; T^2 / 2, T];
Q = q * blkdiag(Theta, Theta);
f = func_calc_air_dens(x_prev, g, beta_tgt);
w = (mvnrnd(zeros(size(F, 1), 1), Q) )';

x_curr = F * x_prev + G * f + G * [0; -g] + w;