function [mse_px, mse_py] = func_BTT_UKF(dyn_param, meas_param, tr_param, beta_UKF, DEMO_FLAG)

%% target
% parameters
g = dyn_param{1};
beta_tgt = dyn_param{2};
q = dyn_param{3};
T = dyn_param{6};
nSteps = 120 / T;
Theta = [T^3 / 3, T^2 / 2; T^2 / 2, T];
Q = q * blkdiag(Theta, Theta);
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
    x_curr = func_BTT_dyn(x_prev, dyn_param);
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
    z = func_rang_bear_meas(pos, meas_param);
    z_r = z(1);
    z_theta = z(2);
    % log
    Z(i, 1) = z_r;
    Z(i, 2) = z_theta;
end

%% local filter
% parameters
sigma_r = meas_param{1};
sigma_theta = meas_param{2};
R = diag([sigma_r^2, sigma_theta^2]);
% initial estimate
x0_UKF = zeros(nStates, 1);
x0_UKF(1) = 230 * 1e3;
x0_UKF(2) = 2300 * cosd(190);
x0_UKF(3) = 90 * 1e3;
x0_UKF(4) = 2300 * sind(190);
P0_UKF = diag([1000^2, 20^2, 1000^2, 20^2]);
x_UKF = zeros(nSteps, nStates);
P_UKF = zeros(nSteps, nStates, nStates);
% unscented Kalman filter
for i = 1 : nSteps
    if i == 1
        x_upd = x0_UKF;
        P_upd = P0_UKF;
    else
        x_upd = (x_UKF(i - 1, :) )';
        P_upd = shiftdim(P_UKF(i - 1, :, :) );
    end
    % predict
    [x_pred, P_pred] = func_UKF_predict(x_upd, P_upd, @func_BTT_dyn, dyn_param, tr_param, Q);
    
    % update
    z = (Z(i, :) )';  
    [x_upd, P_upd] = func_UKF_update(x_pred, P_pred, z, @func_rang_bear_meas, meas_param, tr_param, R);
    
    % log
    x_UKF(i, :) = x_upd';
    P_UKF(i, :, :) = P_upd;
end