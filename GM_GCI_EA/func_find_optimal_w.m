function w_optim = func_find_optimal_w(gm1, gm2, delta)

w_all = 0.01 : delta : 0.99;

J = zeros(length(w_all), 1);

for i = 1 : length(w_all)
    w = w_all(i);
    gm_f = func_GM_GCI_EA(gm1, gm2, w);
    d_1 = func_calc_GM_KLD(gm_f, gm1);
    d_2 = func_calc_GM_KLD(gm_f, gm2);
    J(i) = (d_1 - d_2)^2;
end

[J_min, I] = min(J);

w_optim = w_all(I);