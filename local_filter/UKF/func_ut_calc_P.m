function P = func_ut_calc_P(X)

% read data from sigma-point set
X = sigma.X;
W_m = sigma.W_m;
W_c = sigma.W_c;

L = size(X, 2);

% calculate `P_tfr`
P = [];
for i = 1 : L
    if i == 1
        P = W_c(i) * (sigma(:, i) - m_tfr) * (sigma(:, i) - m_tfr)';
    else
        P = P + W_c(i) * (sigma(:, i) - m_tfr) * (sigma(:, i) - m_tfr)';
    end
end