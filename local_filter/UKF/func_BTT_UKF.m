function [mse_px, mse_py] = func_BTT_UKF(T, q, R, dyn_param, meas_param, ut_param, DEMO_FLAG)

%% target
% read dynamical parameters
g = dyn_param{1};
beta_tgt = dyn_param{2};
F = dyn_param{3};
G = dyn_param{4};
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
    w = (mvnrnd(zeros(nStates, 1), Q) )';
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
    pos_tgt = [x; y];
    z_true = func_rang_bear_meas(pos_tgt, meas_param);
    % measurement noise 
    v = (mvnrnd(zeros(2, 1), R) )';
    % calculate measurement z
    z = z_true + v;
    % log
    Z(i, :) = z';
end

%% local filter
% initial estimate
x0_UKF = zeros(nStates, 1);
x0_UKF(1) = 230 * 1e3;
x0_UKF(2) = 2300 * cosd(190);
x0_UKF(3) = 90 * 1e3;
x0_UKF(4) = 2300 * sind(190);
% x0_UKF = x0_tgt;
P0_UKF = diag([100^2, 10^2, 100^2, 10^2]);
x_UKF = zeros(nSteps, nStates);
P_UKF = zeros(nSteps, nStates, nStates);
% unscented Kalman filter
for i = 1 : nSteps
    i
    if i == 1
        x_upd = x0_UKF;
        P_upd = P0_UKF;
    else
        x_upd = (x_UKF(i - 1, :) )';
        P_upd = shiftdim(P_UKF(i - 1, :, :) );
    end
    % predict
    disp('predict')
    [x_pred, P_pred] = func_UKF_predict(x_upd, P_upd, Q, @func_BTT_dyn, dyn_param, ut_param);
    
    % update
    disp('update')
    z = (Z(i, :) )';  
    [x_upd, P_upd] = func_UKF_update(x_pred, P_pred, z, R, @func_rang_bear_meas, meas_param, ut_param);
    
    % log
    x_UKF(i, :) = x_upd';
    P_UKF(i, :, :) = P_upd;
end

%% calculate MSE
error = x_tgt - x_UKF;
error_px = error(:, 1);
error_py = error(:, 3);
mse_px = error_px .* error_px;
mse_py = error_py .* error_py;

%% DEMO
if DEMO_FLAG == 1
    %% testing
    % target trajectory
    figure
    plot(x_tgt(:, 1), x_tgt(:, 3) )
    title('target trajectory')
    xlabel('x (m)')
    ylabel('y (m)')
    grid on
    % velocity
    figure
    vel = zeros(nSteps, 1);
    for i = 1 : nSteps
        vel(i) = sqrt( (x_tgt(i, 2) )^2 + (x_tgt(i, 4) )^2);
    end
    plot(T * (1 : nSteps), vel);
    xlabel('time (s)')
    ylabel('velocity')
    grid on
    % aerodynamic drag force
    figure
    plot(T * (1 : nSteps), aero_drag_force(:, 1) )
    xlabel('time (s)')
    ylabel('aerodynamic drag force in X (m/s^2)')
    grid on
    figure
    plot(T * (1 : nSteps), aero_drag_force(:, 2) )
    xlabel('time (s)')
    ylabel('aerodynamic drag force in Y (m/s^2)')
    grid on
    figure
    ADF = zeros(nSteps, 1);
    for i = 1 : nSteps
        ADF(i) = sqrt( (aero_drag_force(i, 1) )^2 + (aero_drag_force(i, 2) )^2);
    end
    plot(T * (1 : nSteps), ADF )
    xlabel('time (s)')
    ylabel('aerodynamic drag force (m/s^2)')
    grid on
    axis square
    % range measurement
    figure
    plot(T * (1 : nSteps), Z(:, 1) )
    title('range measurement')
    xlabel('time (s)')
    ylabel('range (m)')
    grid on
    % bearing measurement
    figure
    subplot(2, 1, 1)
    plot(T * (1 : nSteps), Z(:, 2))
    title('bearing measurement')
    xlabel('time (s)')
    ylabel('bearing (rad)')
    grid on
    subplot(2, 1, 2)
    plot(T * (1 : nSteps), rad2deg(Z(:, 2) ) )
    title('bearing measurement')
    xlabel('time (s)')
    ylabel('bearing (deg)')
    grid on
    % unscented Kalman filter
    figure
    plot(x_tgt(:, 1), x_tgt(:, 3) )
    hold on
    plot(x_UKF(:, 1), x_UKF(:, 3), '*')
    legend('true target trajectory', 'UKF')
    hold off
    xlabel('x (m)')
    ylabel('y (m)')
    grid on
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