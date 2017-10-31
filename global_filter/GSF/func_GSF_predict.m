function gm_pred = func_GSF_predict(gm_upd, gm_w, func_handle, dyn_param, ut_param)

N_upd = gm_upd.NumComponents;
N_w = gm_w.NumComponents;

%% Step 1: construct `N_upd` sigma-point sets `sigmas_upd`
for i = 1 : N_upd
    x_upd_i = (gm_upd.mu(i, :) )';
    P_upd_i = gm_upd.Sigma(:, :, i);
    sigma_upd = func_constr_sigma(x_upd_i, P_upd_i, ut_param);
    if i == 1
        sigmas_upd = sigma_upd;
    else
        sigmas_upd = [sigmas_upd; sigma_upd];
    end
end

%% Step 2: unscented transform from `sigma_upd` to `sigma_ut`, construct sigma-point sets `sigmas_ut`
for i = 1 : N_upd
    sigma_tmp = sigmas_upd(i);
    W_m_ut = sigma_tmp.W_m;
    W_c_ut = sigma_tmp.W_c;
    
    X_ut = func_ut_transform(sigma_tmp.X, func_handle, dyn_param);
    
    sigma_ut = struct;
    sigma_ut.X = X_ut;
    sigma_ut.W_m = W_m_ut;
    sigma_ut.W_c = W_c_ut;
    
    if i == 1
        sigmas_ut = sigma_ut;
    else
        sigmas_ut = [sigmas_ut; sigma_ut];
    end
   
end

%% Step 3: calculate `x_pred_i` and `P_pred_ij` based on `sigmas_ut(i)`, construct `x_pred_cell` and `P_pred_cell`
% initial
L = length(sigmas_upd(1).W_m);
x_pred_cell = cell(N_upd, 1);
P_pred_cell = cell(N_upd, N_w);
% calculate `x_pred_i` and `P_pred_ij`
for i = 1 : N_upd
    X_ut = sigmas_ut(i).X;
    W_m = sigmas_ut(i).W_m;
    W_c = sigmas_ut(i).W_c;
    % calculate `x_pred_i`
    for l = 1 : L
        if l == 1
            x_pred_i = W_m(l) * X_ut(:, l);
        else
            x_pred_i = x_pred_i + W_m(l) * X_ut(:, l);
        end
    end
    x_pred_cell{i} = x_pred_i;
    % calculate `P_pred_ij`
    for j = 1 : N_w
        Q_j = gm_w.Sigma(:, :, j);
        for l = 1 : L
            if l == 1
                P_pred_ij = W_c(l) * (X_ut(:, l)- x_pred_i) * (X_ut(:, l)- x_pred_i)' + Q_j;
            else
                P_pred_ij = P_pred_ij + W_c(l) * (X_ut(:, l)- x_pred_i) * (X_ut(:, l)- x_pred_i)' + Q_j;
            end
        end
        P_pred_cell{i, j} = P_pred_ij;
    end  
end

%% Step 4: construct Gaussian mixture model `gm_pred`
% mean
for i = 1 : N_upd
    for j = 1 : N_w
        if i == 1 && j == 1
            x_pred = (x_pred_cell{i})';
        else
            x_pred = [x_pred; (x_pred_cell{i})'];
        end
    end
end
assert(size(x_pred, 1) == N_upd * N_w)
% covariance
for i = 1 : N_upd
    for j = 1 : N_w
        if i == 1 && j == 1
            P_pred = P_pred_cell{i, j};
        else
            P_pred = cat(3, P_pred, P_pred_cell{i, j});
        end
    end
end
assert(size(P_pred, 3) == N_upd * N_w)
% weight
k = 1;
p_pred = zeros(1, N_upd * N_w);
for i = 1 : N_upd
    for j = 1 : N_w
        p_pred(k) = gm_upd.ComponentProportion(i) * gm_w.ComponentProportion(j);
        k = k + 1;
    end
end
assert(abs(sum(p_pred(:) ) - 1) < 1e-2)
% construct Gaussian mixture model
gm_pred = gmdistribution(x_pred, P_pred, p_pred);

%% Step 5: Gaussian mixture reduction (TODO)
