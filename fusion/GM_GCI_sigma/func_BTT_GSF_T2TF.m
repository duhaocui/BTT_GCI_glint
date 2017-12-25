function [mse1_px, mse1_py, mse2_px, mse2_py, mse_f_px, mse_f_py] = func_BTT_GSF_T2TF(T, Q_cell, R_cell, epsilon_w, epsilon_v, dyn_param, meas_param, ut_param, DEMO_FLAG)

%% target
% read dynamical parameters
g = dyn_param{1};
beta_tgt = dyn_param{2};
F = dyn_param{3};
G = dyn_param{4};
% prepare
nSteps = 120 / T;
% initial target state
x0_tgt = [232 * 1e3; 2290 * cosd(190); 88 * 1e3; 2290 * sind(190)];
nStates = length(x0_tgt);
x_tgt = zeros(nSteps, nStates);
aero_drag_force = zeros(nSteps, 2);
% update target state
for i = 1 : nSteps
    if i == 1
        x_prev = x0_tgt;
    else
        x_prev = (x_tgt(i - 1, :) )';
    end
    % process noise (Gaussian mixture model)
    mu_w = zeros(2, size(F, 1) );
    sigma_w = cat(3, Q_cell{1}, Q_cell{2} );
    p_w = [(1 - epsilon_w), epsilon_w];
    gm_w = gmdistribution(mu_w, sigma_w, p_w);
    w = (random(gm_w) )';
    x_curr = func_BTT_dyn(x_prev, dyn_param) + w;
    % log
    x_tgt(i, :) = x_curr';
    f = func_calc_air_dens(x_prev, g, beta_tgt);
    aero_drag_force(i, :) = f';
end

%% measurement
% initial measurement
Z1 = zeros(nSteps, 2);
Z2 = zeros(nSteps, 2);
% update measurement
for i = 1 : nSteps
    x = x_tgt(i, 1);
    y = x_tgt(i, 3);
    pos = [x; y];
    z_true = func_rang_bear_meas(pos, meas_param);   
    % measurement noise (Gaussian mixture model)
    mu_v = zeros(2, 2);
    sigma_v = cat(3, R_cell{1}, R_cell{2});
    p_v = [(1 - epsilon_v), epsilon_v];
    gm_v = gmdistribution(mu_v, sigma_v, p_v);
    v1 = (random(gm_v) )';
    v2 = (random(gm_v) )';
    % calculate measurement z
    z1 = z_true + v1;
    z2 = z_true + v2;
    % log
    Z1(i, :) = z1';
    Z2(i, :) = z2';
end
    
%% Gaussian sum filter (GSF)
% initial filtering estimate (2 Gaussian mixtures, size(F, 1) states)
x0 = zeros(2, size(F, 1) );
x0(1, :) = [230 * 1e3, 2300 * cosd(190), 90 * 1e3, 2300 * sind(190)];
x0(2, :) = [24 * 1e3, 220 * cosd(190), 10 * 1e3, 220 * sind(190)];
P0_cell = cell(2, 1);
P0_cell{1} = diag([1000^2, 20^2, 1000^2, 20^2]);
P0_cell{2} = diag([1500^2, 15^2, 500^2, 25^2]);
P0 = cat(3, P0_cell{1}, P0_cell{2});
gm0 = gmdistribution(x0, P0, [0.6, 0.4]);
% initial filtering log
gm1_cell = cell(nSteps, 1);
gm2_cell = cell(nSteps, 1);
gm_f_cell = cell(nSteps, 1);
% Gaussian sum filtering
for i = 1 : nSteps
    fprintf('current step: %d\n', i)
    if i == 1
        gm1_upd = gm0;
        gm2_upd = gm0;
    else
        gm1_upd = gm1_cell{i - 1};
        gm2_upd = gm2_cell{i - 1};
    end
    % predict
    disp('predict')
    gm1_pred = func_GSF_predict(gm1_upd, gm_w, @func_BTT_dyn, dyn_param, ut_param);
    gm2_pred = func_GSF_predict(gm2_upd, gm_w, @func_BTT_dyn, dyn_param, ut_param);
    % update
    disp('update')
    z1 = (Z1(i, :) )';
    z2 = (Z2(i, :) )';
    [gm1_upd, gm1_upd_full] = func_GSF_update(gm1_pred, z1, gm_v, @func_rang_bear_meas, meas_param, ut_param, i);
    [gm2_upd, gm2_upd_full] = func_GSF_update(gm2_pred, z2, gm_v, @func_rang_bear_meas, meas_param, ut_param, i);
    % T2TF (Chernoff fusion)
    disp('fusion')
    w_Chernoff = 0.5;
    alpha = 1;
    beta = 0;
    kappa = 0;
    ut_param_Chernoff{1} = alpha;
    ut_param_Chernoff{2} = beta;
    ut_param_Chernoff{3} = kappa;
    gm1_power = func_gm_power(gm1_upd_full, w_Chernoff, ut_param_Chernoff);
    gm2_power = func_gm_power(gm2_upd_full, 1 - w_Chernoff, ut_param_Chernoff);
    gm_f = func_Chernoff_gm(gm1_power, gm2_power, w_Chernoff);
    
    % log
    gm1_cell{i} = gm1_upd; 
    gm2_cell{i} = gm2_upd;
    gm_f_cell{i} = gm_f;
    
end

%% calculate MSE
% calculate `x_GSF`
x1_GSF = zeros(nSteps, nStates);
x2_GSF = zeros(nSteps, nStates);
x_f_GSF = zeros(nSteps, nStates);
for i = 1 : nSteps
    gm1 = gm1_cell{i};
    x1_GSF(i, :) = gm1.ComponentProportion * gm1.mu;
    gm2 = gm2_cell{i};
    x2_GSF(i, :) = gm2.ComponentProportion * gm2.mu;
    gm_f = gm_f_cell{i};
    x_f_GSF(i, :) = gm_f.ComponentProportion * gm_f.mu;
end
error1 = x_tgt - x1_GSF;
error1_px = error1(:, 1);
error1_py = error1(:, 3);
mse1_px = error1_px .* error1_px;
mse1_py = error1_py .* error1_py;

error2 = x_tgt - x2_GSF;
error2_px = error2(:, 1);
error2_py = error2(:, 3);
mse2_px = error2_px .* error2_px;
mse2_py = error2_py .* error2_py;

error_f = x_tgt - x_f_GSF;
error_f_px = error_f(:, 1);
error_f_py = error_f(:, 3);
mse_f_px = error_f_px .* error_f_px;
mse_f_py = error_f_py .* error_f_py;

%% plot
if DEMO_FLAG == 1
    % plot estimation error for GSF1
    figure
    plot(T * (1 : nSteps), abs(error1_px) )
    xlabel('time (s)')
    ylabel('estimation error in X (absolute value)')
    grid on
    figure
    plot(T * (1 : nSteps), abs(error1_py) )
    ylabel('estimation error in Y (absolute value)')
    xlabel('time (s)')
    grid on
    % plot estimation error for GSF2
    figure
    plot(T * (1 : nSteps), abs(error2_px) )
    xlabel('time (s)')
    ylabel('estimation error in X (absolute value)')
    grid on
    figure
    plot(T * (1 : nSteps), abs(error2_py) )
    ylabel('estimation error in Y (absolute value)')
    xlabel('time (s)')
    grid on
    % plot estimation error for fusion
    figure
    plot(T * (1 : nSteps), abs(error_f_px) )
    xlabel('time (s)')
    ylabel('estimation error in X (absolute value)')
    grid on
    figure
    plot(T * (1 : nSteps), abs(error_f_py) )
    ylabel('estimation error in Y (absolute value)')
    xlabel('time (s)')
    grid on
end