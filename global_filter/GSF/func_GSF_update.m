function gm_upd = func_GSF_update(gm_pred, z, gm_v, func_handle, meas_param, ut_param)

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
    
    X_ut = func_ut_transform(sigma_tmp.X, func_handle, meas_param);
    
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

%% Step 3: calculate `z_pred_i`, `P_xz_pred_i` and `P_z_pred_ij` based on `sigmas_ut(i)`, `sigma_pred(i)`, construct `z_pred_cell`, `P_xz_pred_cell` and `P_z_pred_cell`
% initial
L = length(sigmas_pred(1).W_m);
z_pred_cell = cell(N_pred, 1);
P_xz_pred_cell = cell(N_pred, 1);
P_z_pred_cell = cell(N_pred, N_v);
% calculate `z_pred_i`, `P_xz_pred_i` and `P_z_pred_ij`
for i = 1 : N_pred
    Z_ut = sigmas_ut(i).X;
    X_pred = sigmas_pred(i).X;
    W_m = sigmas_ut(i).W_m;
    W_c = sigmas_ut(i).W_c;
    x_pred_i = (gm_pred.mu(i, :) )';
    % calculate `z_pred_i`
    for l = 1 : L
        if l == 1
            z_pred_i = W_m(l) * Z_ut(:, l);
        else
            z_pred_i = z_pred_i + W_m(l) * Z_ut(:, l);
        end
    end
    z_pred_cell{i} = z_pred_i;
    % calculate `P_xz_pred_i`
    for l = 1 : L
        if l == 1
            P_xz_pred_i = W_c(l) * (X_pred(:, l) - x_pred_i) * (Z_ut(:, l) - z_pred_i)'; 
        else
            P_xz_pred_i = P_xz_pred_i + W_c(l) * (X_pred(:, l) - x_pred_i) * (Z_ut(:, l) - z_pred_i)'; 
        end
    end
    P_xz_pred_cell{i} = P_xz_pred_i;
    % calculate `P_z_pred_ij`
    for j = 1 : N_v
        R_j = gm_v.Sigma(:, :, j);
        for l = 1 : L
            if l == 1
                P_z_pred_ij = W_c(l) * (Z_ut(:, l) - z_pred_i) * (Z_ut(:, l) - z_pred_i)' + R_j;
            else
                P_z_pred_ij = P_z_pred_ij + W_c(l) * (Z_ut(:, l) - z_pred_i) * (Z_ut(:, l) - z_pred_i)' + R_j;
            end
        end
        P_z_pred_cell{i, j} = P_z_pred_ij;
    end
    
end


%% Step 4:

%% Step 5: