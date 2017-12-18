function gm = func_constr_gm(mean, cov, w)

gm = struct;

gm.x = mean;
gm.P = cov;
gm.w = w;

assert(sum(w) - 1 < 1e-1)