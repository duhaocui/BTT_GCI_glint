function [X, W] = func_constr_sigma(m, P, tr_param)
% ref. [1]: Gaussian sum unscented kalman filter with adaptive scaling parameters
% ref. [2]: Radar data processing with applications (He You), P72

n_x = length(m);
kappa = tr_param{1};

% construct sigma points
X = zeros(n_x, 2 * n_x + 1);
X(:, 1) = m;

% P_sqrt = chol(P)';
P_sqrt = schol(P)';

for i = 1 : n_x
    X(:, i + 1) = m + sqrt(n_x + kappa) * P_sqrt(:, i);
    X(:, i + 1 + n_x) = m - sqrt(n_x + kappa) * P_sqrt(:, i);
end

% calculate `W`
W = zeros(2 * n_x + 1, 1);
W(1) = kappa / (n_x + kappa);
for i = 1 : n_x
    W(i + 1) = 1 / (2 * (n_x + kappa) );
end

for i = 1 : n_x
    W(i + 1 + n_x) = W(i + 1);
end


