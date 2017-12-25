function gm_f = func_Chernoff_gm(gm1_power, gm2_power, w)

gm1_power_ea = func_tsfr_gm_ea(gm1_power);
gm2_power_ea = func_tsfr_gm_ea(gm2_power);
gm_f_EA = func_GM_GCI_EA(gm1_power_ea, gm2_power_ea, w);

%% construct Gaussian mixture model
% initial
k = length(gm_f_EA.w);
d = length(gm_f_EA.x{1});
mu = zeros(k, d);
sigma = zeros(d, d, k);
p = zeros(1, k);
% construct GMM
% get `mu`
for i = 1 : k
    mu(i, :) = (gm_f_EA.x{i})';
end
% get `sigma`
for i = 1 : k
   sigma(:, :, k) = gm_f_EA.P{i};
end
% get `p`
for i = 1 : k
    p(i) = gm_f_EA.w(i);
end

gm_f = gmdistribution(mu, sigma, p);