function gm_upd = func_GSF_update(gm_pred, z, gm_v, func_handle, meas_param, ut_param, step_curr)

N_pred = gm_pred.NumComponents;
N_v = gm_v.NumComponents;

%% Step 1: construct `N_pred` sigma-point sets `sigmas_pred`
for i = 1 : N_pred
    x_pred_i = (gm_pred.mu(i, :) )';
    P_pred_i = gm_pred.Sigma(:, :, i);
    sigma_pred = func_constr_sigma(x_pred_i, P_pred_i, ut_param);
    if i == 1
        sigmas_pred = sigma_pred;
    else
        sigmas_pred = [sigmas_pred; sigma_pred];
    end
end

%% Step 2: unscented transform from `sigma_pred` to `sigma_ut`, construct sigma-point sets `sigmas_ut`
for i = 1 : N_pred
    sigma_tmp = sigmas_pred(i);
    W_m_ut = sigma_tmp.W_m;
    W_c_ut = sigma_tmp.W_c;
    
    X_ut_i = func_ut_transform(sigma_tmp.X, func_handle, meas_param);
    
    sigma_ut = struct;
    sigma_ut.X = X_ut_i;
    sigma_ut.W_m = W_m_ut;
    sigma_ut.W_c = W_c_ut;
    
    if i == 1
        sigmas_ut = sigma_ut;
    else
        sigmas_ut = [sigmas_ut; sigma_ut];
    end
end

%% Step 3: calculate `z_pred_i`, `P_xz_pred_i`, `P_z_pred_ij` and `phi_z`
% based on `sigmas_ut(i)`, `sigma_pred(i)`, construct `z_pred_cell`, `P_xz_pred_cell` and `P_z_pred_cell`
% initial
L = length(sigmas_pred(1).W_m);
z_pred_cell = cell(N_pred, 1);
P_xz_pred_cell = cell(N_pred, 1);
P_z_pred_cell = cell(N_pred, N_v);
phi_z = zeros(N_pred, N_v);
% calculate `z_pred_i`, `P_xz_pred_i`, `P_z_pred_ij` and `phi_z`
for i = 1 : N_pred
    Z_ut_i = sigmas_ut(i).X;
    X_pred_i = sigmas_pred(i).X;
    W_m = sigmas_ut(i).W_m;
    W_c = sigmas_ut(i).W_c;
    x_pred_i = (gm_pred.mu(i, :) )';
    % calculate `z_pred_i`
    for l = 1 : L
        if l == 1
            z_pred_i = W_m(l) * Z_ut_i(:, l);
        else
            z_pred_i = z_pred_i + W_m(l) * Z_ut_i(:, l);
        end
    end
    z_pred_cell{i} = z_pred_i;
    % calculate `P_xz_pred_i`
    for l = 1 : L
        if l == 1
            P_xz_pred_i = W_c(l) * (X_pred_i(:, l) - x_pred_i) * (Z_ut_i(:, l) - z_pred_i)'; 
        else
            P_xz_pred_i = P_xz_pred_i + W_c(l) * (X_pred_i(:, l) - x_pred_i) * (Z_ut_i(:, l) - z_pred_i)'; 
        end
    end
    P_xz_pred_cell{i} = P_xz_pred_i;
    % calculate `P_z_pred_ij`
    for j = 1 : N_v
        R_j = gm_v.Sigma(:, :, j);
        for l = 1 : L
            if l == 1
                P_z_pred_ij = W_c(l) * (Z_ut_i(:, l) - z_pred_i) * (Z_ut_i(:, l) - z_pred_i)';
            else
                P_z_pred_ij = P_z_pred_ij + W_c(l) * (Z_ut_i(:, l) - z_pred_i) * (Z_ut_i(:, l) - z_pred_i)';
            end
        end
        P_z_pred_ij = P_z_pred_ij + R_j;
        P_z_pred_cell{i, j} = P_z_pred_ij;
    end
    % calculate `phi_z`
    for j = 1 : N_v
        phi_z(i, j) = mvnpdf(z, z_pred_cell{i}, P_z_pred_cell{i, j});
    end
end

%% Step 4: unscented Kalman filtering (update phase)
x_upd_cell = cell(N_pred, N_v);
P_upd_cell = cell(N_pred, N_v);
for i = 1 : N_pred
    x_pred_i = (gm_pred.mu(i, :) )';
    z_pred_i = z_pred_cell{i};
    P_pred_i = gm_pred.Sigma(:, :, i);
    for j = 1 : N_v
        K = P_xz_pred_cell{i} / P_z_pred_cell{i, j};
        x_upd_ij = x_pred_i + K * (z - z_pred_i);
        P_upd_ij = P_pred_i - K * P_z_pred_cell{i, j} * K';
        % log
        x_upd_cell{i, j} = x_upd_ij;
        P_upd_cell{i, j} = P_upd_ij;
    end
end

%% Step 5: construct Gaussian mixture model `gm_upd`
% mean
for i = 1 : N_pred
    for j = 1 : N_v
        if i == 1 && j == 1
            x_upd = (x_upd_cell{i, j})';
        else
            x_upd = [x_upd; (x_upd_cell{i, j})'];
        end
    end
end
assert(size(x_upd, 1) == N_pred * N_v)
% covariance
for i = 1 : N_pred
    for j = 1 : N_v
        if i == 1 && j == 1
            P_upd = P_upd_cell{i, j};
        else
            P_upd = cat(3, P_upd, P_upd_cell{i, j});
        end
    end
end
assert(size(P_upd, 3) == N_pred * N_v)
% weight
k = 1;
p_upd = zeros(1, N_pred * N_v);
for i = 1 : N_pred
    for j = 1 : N_v
        p_upd(k) = gm_pred.ComponentProportion(i) * gm_v.ComponentProportion(j) * phi_z(i, j);
        k = k + 1;
    end
end
p_upd = p_upd / sum(p_upd);
assert(abs(sum(p_upd(:) ) - 1) < 1e-2 )
% prune
elim_threshold = 1e-5;
w_old = p_upd';
x_old = x_upd';
P_old = P_upd;
[w_new, x_new, P_new]= gaus_prune(w_old, x_old, P_old, elim_threshold);
w_new = w_new / sum(w_new);
assert(abs(sum(w_new) - 1) < 1e-2)
% construct Gaussian mixture model
gm_upd = gmdistribution(x_new', P_new, w_new');

if step_curr == 40
    save('gm_upd.mat', 'gm_upd');
end

%% Step 6: Gaussian mixture management
% merge
gm_upd = gm_merge(gm_upd);
% cap
gm_upd = gm_cap(gm_upd);
