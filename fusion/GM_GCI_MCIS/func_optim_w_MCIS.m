function w_optim = func_optim_w_MCIS(gm1, gm2, delta)

w_all = 0.01 : delta : 0.99;
J = zeros(length(w_all), 1);

for i = 1 : length(w_all)
    w = w_all(i);
    gm_f = func_Chernoff_approx(gm1, gm2, w);

    gm1_ea = func_tsfr_gm_ea(gm1);
    gm2_ea = func_tsfr_gm_ea(gm2);
    gm_f_ea = func_tsfr_gm_ea(gm_f);
    
    d_1 = func_calc_GM_KLD(gm_f_ea, gm1_ea);
    d_2 = func_calc_GM_KLD(gm_f_ea, gm2_ea);
    
    J(i) = (d_1 - d_2)^2;
end

[J_min, I] = min(J);
w_optim = w_all(I);