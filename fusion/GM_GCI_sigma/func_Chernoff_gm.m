function gm_f = func_Chernoff_gm(gm1_power, gm2_power, w)

gm1_power_ea = func_tsfr_gm_ea(gm1_power);
gm2_power_ea = func_tsfr_gm_ea(gm2_power);
gm_f = func_GM_GCI_EA(gm1_power_ea, gm2_power_ea, w);