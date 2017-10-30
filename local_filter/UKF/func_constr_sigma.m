function [X, W_m, W_c] = func_constr_sigma(m, P, tr_param)
% ref. [1]: Gaussian sum unscented kalman filter with adaptive scaling parameters
% ref. [2]: Radar data processing with applications (He You), P72

n_x = length(m);
alpha = tr_param{1};
beta = tr_param{2};
kappa = tr_param{3};

% calculate lambda
lambda = alpha^2 * (n_x + kappa) - n_x;

%% calculate sigma points
% initial
X = zeros(n_x, 2 * n_x + 1);
% calcualte
tmp = chol( (n_x + lambda) * P, 'lower');

% testing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tmp = schol( (n_x + lambda) * P);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X(:, 1) = m;
for i = 1 : n_x
    X(:, i + 1) = m + tmp(:, i);
    X(:, i + 1 + n_x) = m - tmp(:, i);
end

%% calcualte W_m
% initial
W_m = zeros(2 * n_x + 1, 1);
% calculate
W_m(1) = lambda / (n_x + lambda);
W_m(2 : end) = 1 / (2 * (n_x + lambda) );

% calculate W_c
% initial
W_c = zeros(2 * n_x + 1, 1);
% calculate
W_c(1) = lambda / (n_x + lambda) + (1 - alpha^2 + beta);
W_c(2 : end) = 1 / (2 * (n_x + lambda) );