function gm_ea = func_tsfr_gm_ea(gm)

mean_cell = cell(gm.NumComponents, 1);
cov_cell = cell(gm.NumComponents, 1);
w = zeros(gm.NumComponents, 1);

for i = 1:gm.NumComponents
    mean_cell{i} = (gm.mu(i, :) )';
    cov_cell{i} = gm.Sigma(:, :, i);
    w(i) = gm.ComponentProportion(i);
end

gm_ea = func_constr_gm(mean_cell, cov_cell, w);