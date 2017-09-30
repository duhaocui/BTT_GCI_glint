function [X, W_m, W_c] = func_constr_sigma(m, P, tr_param)
% ref. [1]: Optimal Filterign with Kalman Filters and Smoothers, P29
% ref. [2]: ??????????

n = length(m);
alpha = tr_param{1};
beta = tr_param{2};
kappa = tr_param{3};

% calculate `lambda`
lambda = alpha^2 * (n + kappa) - n;

% construct sigma points
X = zeros(n, 2 * n + 1);
X(:, 1) = m;

% P_sqrt = chol(P)';
P_sqrt = schol(P)';

for i = 1 : n
    X(:, i + 1) = m + sqrt(n + lambda) * P_sqrt(:, i);
    X(:, i + 1 + n) = m - sqrt(n + lambda) * P_sqrt(:, i);
end

% calculate `W_m`
W_m = zeros(2 * n + 1, 1);
W_m(1) = lambda / (n + lambda);
for i = 1 : 2 * n
    W_m(1 + i) = 1 / (2 * (n + lambda) );
end

% calculate `W_c`
W_c = zeros(2 * n + 1, 1);
W_c(1) = lambda / (n + lambda) + (1 - alpha^2 + beta);
for i = 1 : 2 * n
    W_c(1 + i) = 1 / (2 * (n + lambda) );
end


