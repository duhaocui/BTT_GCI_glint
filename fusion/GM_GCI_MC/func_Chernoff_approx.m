function gm_f = func_Chernoff_approx(gm1, gm2, w)

N_a = gm1.NumComponents;
N_b = gm2.NumComponents;

p = gm1.ComponentProportion;
q = gm2.ComponentProportion;

a = gm1.mu;
b = gm2.mu;

A = gm1.Sigma;
B = gm2.Sigma;

C_cell = cell(N_a, N_b);
c_cell = cell(N_a, N_b);
r = zeros(N_a, N_b);

%% calculate `den_r` for `r_ij`
for k = 1 : N_a
    p_k = p(k);
    for l = 1 : N_b
        q_l = q(l);
        if k == 1 && l == 1
            den_r = (p_k^w) * (q_l^(1 - w) );
        else
            den_r = den_r + (p_k^w) * (q_l^(1 - w) );
        end
    end
end

%% calculate `C_ij`, `c_ij`, `r_ij`
for i = 1 : N_a
    A_i = A(:, :, i);
    a_i = (a(i, :) )';
    p_i = p(i);
    for j = 1 : N_b
        B_j = B(:, :, j);
        b_j = (b(j, :) )';
        q_j = q(j);
        inv_C_ij = w * inv(A_i) + (1 - w) * inv(B_j);
        C_ij = inv(inv_C_ij);
        c_ij = C_ij * (w * inv(A_i) * a_i + (1 - w) * inv(B_j) * b_j);        
        r_ij = (p_i^w) * (q_j^(1 - w) ) / den_r;
        % log
        C_cell{i, j} = C_ij;
        c_cell{i, j} = c_ij;
        r(i, j) = r_ij;
    end
end

%% construct Gaussian mixture
% initial
k = N_a * N_b;
d = size(a, 2);
mu = zeros(k, d);
sigma = zeros(d, d, k);
p = zeros(1, k);
% construct GM
k = 1;
for i = 1 : N_a
    for j = 1 : N_b
        mu(k, :) = (c_cell{i, j})';
        sigma(:, :, k) = C_cell{i, j};
        p(k) = r(i, j);
        k = k + 1;
    end
end
gm_f = gmdistribution(mu, sigma, p);
