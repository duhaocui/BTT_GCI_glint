function [d, klfg] = func_calc_GM_KLD(gm1, gm2)

% Outputs:
%   d             the approximate KL divergence D(f||g)
%   klfg(kf,kg)   the exact KL divergence between the components of f and g

kf = length(gm1.x);
kg = length(gm2.x);

assert(length(gm1.x{1} ) == length(gm2.x{1} ))
p = length(gm1.x{1} );

mf = zeros(kf, p);
vf = zeros(p, p, kf);
wf = zeros(kf, 1);

mg = zeros(kg, p);
vg = zeros(p, p, kg);
wg = zeros(kg, 1);

for i = 1 : kf
    mf(i, :) = (gm1.x{i} )';
    vf(:, :, i) = gm1.P{i};
    wf(i) = gm1.w(i);
end

for i = 1 : kg
    mg(i, :) = (gm2.x{i} )';
    vg(:, :, i) = gm2.P{i};
    wg(i) = gm2.w(i);
end

[d, klfg] = gaussmixk(mf, vf, wf, mg, vg, wg);