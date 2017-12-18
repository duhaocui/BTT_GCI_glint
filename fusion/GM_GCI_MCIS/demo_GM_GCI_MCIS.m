clear
clc
close all

load('gm1.mat') 
gm1 = gm_upd;
load('gm2.mat') 
gm2 = gm_upd;

% first order approximation of Chernoff fusion
delta = 0.1;
w_optim = func_optim_w_MCIS(gm1, gm2, delta);
gm_f = func_Chernoff_approx(gm1, gm2, w_optim);