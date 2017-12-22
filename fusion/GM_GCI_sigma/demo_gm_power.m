clear
clc
close all

load('gm1.mat') 
gm1 = gm_upd;

w = 0.5;

alpha = 1;
beta = 0;
kappa = 0;
ut_param{1} = alpha;
ut_param{2} = beta;
ut_param{3} = kappa;

gm_power = func_gm_power(gm1, w, ut_param);