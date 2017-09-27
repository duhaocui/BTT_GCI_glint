function [x, P] = func_GM2Gaussian(gm)

N = length(gm.x);

% calculatate m 
for i = 1 : N
    if i == 1
        x = gm.w(i) * gm.x{i};
    else
        x = x + gm.w(i) * gm.x{i};
    end
end

% calculate P
for i = 1 : N
    if i == 1
        P = gm.w(i) * (gm.P{i} + (gm.x{i} - x) * (gm.x{i} - x)' );
    else
        P = P + gm.w(i) * (gm.P{i} + (gm.x{i} - x) * (gm.x{i} - x)' );
    end
end

