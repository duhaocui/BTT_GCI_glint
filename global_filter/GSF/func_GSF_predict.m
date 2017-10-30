function gm_pred = func_GSF_predict(gm_upd, func_handle, dyn_param, ut_param)

N_k_k = gm_upd.NumComponents;

%% Step 1: construct `N_k_k` sigma-point sets `sigmas`
for i = 1 : N_k_k
    m = (gm_upd.mu(i, :) )';
    P = gm_upd.Sigma(:, :, i);
    sigma_upd = func_constr_sigma(m, P, ut_param);
    if i == 1
        sigmas_upd = sigma_upd;
    else
        sigmas_upd = [sigmas_upd; sigma_upd];
    end
end

%% Step 2: unscented from `sigma_upd` to `sigma_pred`
t = 1;