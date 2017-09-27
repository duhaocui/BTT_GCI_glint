clear
clc 
close all

filename = 'gm_data.mat';
myVars = {'gm1', 'gm2'};
load(filename, myVars{:})

[d, klfg] = func_calc_GM_KLD(gm1, gm2);