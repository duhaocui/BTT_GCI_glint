function gm_new = gm_merge(gm_old)

merge_threshold = 4;  

w_old = (gm_old.ComponentProportion)';
x_old = (gm_old.mu)';
P_old = gm_old.Sigma;

[w_new, x_new, P_new]= gaus_merge(w_old, x_old, P_old, merge_threshold);
w_new = w_new / sum(w_new);
assert(abs(sum(w_new) - 1) < 1e-2)

gm_new = gmdistribution(x_new', P_new, w_new);