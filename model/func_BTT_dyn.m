function x_curr = func_BTT_dyn(x_prev, dyn_param)

g = dyn_param{1};
beta_tgt = dyn_param{2};
F = dyn_param{3};
G = dyn_param{4};

f = func_calc_air_dens(x_prev, g, beta_tgt);

x_curr = F * x_prev + G * f + G * [0; -g];