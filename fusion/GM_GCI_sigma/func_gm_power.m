function gm_power = func_gm_power(gm, w, ut_param)

%% initial
M = gm.NumComponents;
n_x = gm.NumVariables;
L = 2 * n_x + 1;
alpha = zeros(M, 1);
pi = zeros(M, L);
X_sigma_cell = cell(M, 1);

%% get `alpha`
for i = 1 : M
    alpha(i) = gm.ComponentProportion(i);
end

%% get `pi` and `X_cell`
for i = 1 : M
    mu = (gm.mu(i, :) )';
    P = gm.Sigma(:, :, i);
    % construct sigma points
    sigma = func_constr_sigma(mu, P, ut_param);
    W_m = sigma.W_m;
    assert(length(W_m) == L)
    % get `pi`
    for l = 1 : L
        pi(i, l) = W_m(l);
    end
    % get `X_sigma_cell`
    X_sigma_cell{i} = sigma.X;  
end

%% initial
A = zeros(M * L, M);
b = zeros(M * L, 1);
W = zeros(M * L, M * L);

%% get `A`
for i = 1 : M
    for l = 1 : L
        for m = 1 : M
            x_sigma_i_l = (X_sigma_cell{i}(:, l) )';
            x_m = gm.mu(m, :);
            P_m = gm.Sigma(:, :, m);
            A(L * (i - 1) + l, m) = mvnpdf(x_sigma_i_l, x_m, P_m / w);
        end
    end
end

%% get `b`
for i = 1 : M
    for l = 1 : L
        x = (X_sigma_cell{i}(:, l) )';
        b(L * (i - 1) + l) = (pdf(gm, x) )^w;
    end
end

%% get `W`
for i = 1 : M
    for l = 1 : L
        W(L * (i - 1) + l, L * (i - 1) + l) = alpha(i) * pi(i, l);
    end
end

%% solve nonnegative linear least-squares problem with `lsqnonneg()`
C = sqrt(W) * A;
d = sqrt(W) * b;
beta = lsqnonneg(C, d);

%% construct Gaussian mixture model
% get `mu_power`
mu_power = gm.mu;
% get `sigma_power`
sigma_power = gm.Sigma / w;
% get `p_power`
p_power = zeros(gm.NumComponents, 1);
for i = 1 : gm.NumComponents
    p_power(i) = beta(i) / sum(beta);
end
% prune
elim_threshold = 1e-5;
p_power_old = p_power;
mu_power_old = mu_power';
sigma_power_old = sigma_power;
[p_power_new, mu_power_new, sigma_power_new]= gaus_prune(p_power_old, mu_power_old, sigma_power_old, elim_threshold);
p_power_new = p_power_new / sum(p_power_new);
assert(abs(sum(p_power_new) - 1) < 1e-1)
% construct Gaussian mixture model
gm_power = gmdistribution(mu_power_new', sigma_power_new, p_power_new');