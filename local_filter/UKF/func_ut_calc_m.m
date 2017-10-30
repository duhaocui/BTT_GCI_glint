function m = func_ut_calc_m(X, W_m, W_c)

L = size(X, 2);

% calculate `m_tfr`
m = [];
for i = 1 : L
    if i == 1
        m = W_m(i) * X(:, i);
    else
        m = m + W_m(i) * X(:, i);
    end
end