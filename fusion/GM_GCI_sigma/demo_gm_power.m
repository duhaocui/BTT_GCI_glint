clear
clc
close all

load('gm1.mat') 
gm1 = gm_upd;
load('gm2.mat')
gm2 = gm_upd;

w = 0.5;

alpha = 1;
beta = 0;
kappa = 0;
ut_param{1} = alpha;
ut_param{2} = beta;
ut_param{3} = kappa;

gm1_power = func_gm_power(gm1, w, ut_param);
gm2_power = func_gm_power(gm2, 1 - w, ut_param);

gm_f = func_Chernoff_gm(gm1_power, gm2_power, w);