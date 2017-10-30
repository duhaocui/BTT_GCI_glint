function [mse_px, mse_py] = func_BTT_EKF(T, q, g, beta_tgt, F, G, x_R, y_R, sigma_r, sigma_theta, beta_EKF, DEMO_FLAG)

%% target
% parameters
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
    f = func_calc_air_dens(x_prev, g, beta_tgt);
    w = (mvnrnd(zeros(nStates, 1), Q) )';
    x_curr = F * x_prev + G * f + G * [0; -g] + w;
    % log
    x_tgt(i, :) = x_curr';
    aero_drag_force(i, :) = f';
end

%% measurement
% initial measurement
Z = zeros(nSteps, 2);
% update measurement
for i = 1 : nSteps
    x = x_tgt(i, 1);
    y = x_tgt(i, 3);
    v_r = normrnd(0, sigma_r);
    v_theta = normrnd(0, sigma_theta); 
    z_r = sqrt( (x - x_R)^2 + (y - y_R)^2) + v_r;
    z_theta = atan2( (y - y_R), (x - x_R) ) + v_theta;
    % log
    Z(i, 1) = z_r;
    Z(i, 2) = z_theta;
end

%% local filter
% parameters
R = diag([sigma_r^2, sigma_theta^2]);
% initial estimate
x0_EKF = zeros(nStates, 1);
x0_EKF(1) = 230 * 1e3;
x0_EKF(2) = 2300 * cosd(190);
x0_EKF(3) = 90 * 1e3;
x0_EKF(4) = 2300 * sind(190);
P0_EKF = diag([1000^2, 20^2, 1000^2, 20^2]);
x_EKF = zeros(nSteps, nStates);
P_EKF = zeros(nSteps, nStates, nStates);
% extended Kalman filter
for i = 1 : nSteps
    if i == 1
        x_upd = x0_EKF;
        P_upd = P0_EKF;
    else
        x_upd = (x_EKF(i - 1, :) )';
        P_upd = shiftdim(P_EKF(i - 1, :, :) );
    end
    % predict
    f_J = func_calc_f_J(g, beta_EKF, x_upd);
    [x_pred, P_pred] = func_EKF_pred(x_upd, P_upd, F, G, g, f_J, Q, beta_EKF);
    % update
    z = (Z(i, :) )';
    h_J = func_calc_h_J(x_R, y_R, x_pred);
    [x_upd, P_upd] = func_EKF_upd(x_R, y_R, x_pred, P_pred, z, h_J, R);
    % log
    x_EKF(i, :) = x_upd';
    P_EKF(i, :, :) = P_upd;
end

%% calculate MSE
error = x_tgt - x_EKF;
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
    % extended Kalman filter
    figure
    plot(x_tgt(:, 1), x_tgt(:, 3) )
    hold on
    plot(x_EKF(:, 1), x_EKF(:, 3), '*')
    legend('true target trajectory', 'EKF')
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