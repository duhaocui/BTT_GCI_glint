function [mse_px, mse_py] = func_BTT_GSF(T, Q_cell, R_cell, epsilon_w, epsilon_v, dyn_param, meas_param, ut_param, DEMO_FLAG)

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
Z = zeros(nSteps, 2);
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
    v = (random(gm_v) )';
    % calculate measurement z
    z = z_true + v;
    % log
    Z(i, :) = z';
end
    
%% Gaussian sum filter
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
gm_cell = cell(nSteps, 1);
% Gaussian sum filtering
for i = 1 : nSteps
    fprintf('current step: %d\n', i)
    if i == 1
        gm_upd = gm0;
    else
        gm_upd = gm_cell{i - 1};
    end
    % predict
    disp('predict')
    gm_pred = func_GSF_predict(gm_upd, gm_w, @func_BTT_dyn, dyn_param, ut_param);
    % update
    disp('update')
    z = (Z(i, :) )';
    gm_upd = func_GSF_update(gm_pred, z, gm_v, @func_rang_bear_meas, meas_param, ut_param, i);
    % log
    gm_cell{i} = gm_upd; 
end

%% calculate MSE
% calculate `x_GSF`
x_GSF = zeros(nSteps, nStates);
for i = 1 : nSteps
    gm = gm_cell{i};
    x_GSF(i, :) = gm.ComponentProportion * gm.mu;
end
error = x_tgt - x_GSF;
error_px = error(:, 1);
error_py = error(:, 3);
mse_px = error_px .* error_px;
mse_py = error_py .* error_py;

%% plot
if DEMO_FLAG == 1
    % plot estimation error
    figure
    plot(T * (1 : nSteps), abs(error_px) )
    xlabel('time (s)')
    ylabel('estimation error in X (absolute value)')
    grid on
    figure
    plot(T * (1 : nSteps), abs(error_py) )
    ylabel('estimation error in Y (absolute value)')
    xlabel('time (s)')
    grid on
end