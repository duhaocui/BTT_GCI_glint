function [mse_px, mse_py, beta_est] = func_BTT_IMM(T, q, g, beta_tgt, F, G, x_R, y_R, sigma_r, sigma_theta, nModes, beta_IMM, DEMO_FLAG)

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
p = zeros(nModes, nModes);
if nModes == 2
    p(1, 1) = 0.9;
    p(1, 2) = 1 - p(1, 1);
    p(2, 1) = p(1, 2);
    p(2, 2) = p(1, 1);
elseif nModes == 3
    p = [0.9, 0.05, 0.05;
         0.05, 0.9, 0.05;
         0.05, 0.05, 0.9];
end
% initial estimate
x0_IMM = zeros(nStates, 1);
x0_IMM(1) = 230 * 1e3;
x0_IMM(2) = 2300 * cosd(190);
x0_IMM(3) = 90 * 1e3;
x0_IMM(4) = 2300 * sind(190);
P0_IMM = diag([1000^2, 20^2, 1000^2, 20^2]);
mu0 = ones(nModes, 1) / nModes;
x_IMM_prior = cell(nSteps, nModes);
P_IMM_prior = cell(nSteps, nModes);
x_IMM = cell(nSteps, nModes);
P_IMM = cell(nSteps, nModes);
mu = zeros(nSteps, nModes);
% IMM
for i = 1 : nSteps
    % step 1
    if i == 1
        mu_km1 = mu0;
    else
        mu_km1 = (mu(i - 1, :) )';
    end
    % calculate the normalization constants c_bar
    c_bar = zeros(nModes, 1);
    for j = 1 : nModes
        for k = 1 : nModes
            if k == 1
                c_bar(j) = p(k, j) * mu_km1(k);
            else
                c_bar(j) = c_bar(j) + p(k, j) * mu_km1(k);
            end
        end
    end
    % calculate the mixing probabilities
    mu_cond_km1 = zeros(nModes, nModes);
    for j = 1 : nModes
        for k = 1 : nModes
            mu_cond_km1(j, k) = p(j, k) * mu_km1(j) / c_bar(k);
        end
    end
    
    % step 2
    % get previous estimate matched to mode j
    x_km1 = cell(nModes, 1);
    P_km1 = cell(nModes, 1);
    for j = 1 : nModes
        if i == 1
            x_km1{j} = x0_IMM;
            P_km1{j} = P0_IMM;
        else
            x_km1{j} = x_IMM{i - 1, j};
            P_km1{j} = P_IMM{i - 1, j};
        end
    end
    % initial mixed initial condition
    x0_km1 = cell(nModes, 1);
    for j = 1 : nModes
        for k = 1 : nModes
            if k == 1
                x0_km1{j} = mu_cond_km1(k, j) * x_km1{k};
            else
                x0_km1{j} = x0_km1{j} + mu_cond_km1(k, j) * x_km1{k};
            end
        end
    end
    P0_km1 = cell(nModes, 1);
    for j = 1 : nModes
        for k = 1 : nModes
            if k == 1
                P0_km1{j} = mu_cond_km1(k, j) * (P_km1{k} + (x_km1{k} - x0_km1{j} ) * (x_km1{k} - x0_km1{j} )' );
            else
                P0_km1{j} = P0_km1{j} + mu_cond_km1(k, j) * (P_km1{k} + (x_km1{k} - x0_km1{j} ) * (x_km1{k} - x0_km1{j} )' );
            end
        end
    end
    
    
    % step 3
    % mode-matched filtering
    z = (Z(i, :) )';
    Lambda = zeros(nModes, 1);
    for j = 1 : nModes
        % predict
        f_J = func_calc_f_J(g, beta_IMM(j), x0_km1{j});
        [x_pred, P_pred] = func_EKF_pred(x0_km1{j}, P0_km1{j}, F, G, g, f_J, Q, beta_IMM(j) );
        % update
        h_J = func_calc_h_J(x_R, y_R, x_pred);
        [x_upd, P_upd] = func_EKF_upd(x_R, y_R, x_pred, P_pred, z, h_J, R);
        % calculate the likelihood functions matched to mode j
        z_r = sqrt( (x_pred(1) - x_R)^2 + (x_pred(3) - y_R)^2);
        z_theta = atan2( (x_pred(3) - y_R), (x_pred(1) - x_R) );
        mean_prior = [z_r; z_theta];
        cov_prior = h_J * P_pred * h_J' + R; 
        Lambda(j) = mvnpdf(z', mean_prior', cov_prior); 
        % log
        x_IMM_prior{i, j} = x_pred;
        P_IMM_prior{i, j} = P_pred;
        x_IMM{i, j} = x_upd;
        P_IMM{i, j} = P_upd;   
    end
    
    % step 4
    % calculate the normalization constants c
    for j = 1 : nModes
        if j == 1
            c = Lambda(j) * c_bar(j);
        else
            c = c + Lambda(j) * c_bar(j);
        end
    end
    % mode probability update
    for j = 1 : nModes
        mu(i, j) = Lambda(j) * c_bar(j) / c;
    end
    
end

% step 5
x_comb_IMM = zeros(nSteps, nStates);
P_comb_IMM = zeros(nSteps, nStates, nStates);
% estimate and covariance combination
for i = 1 : nSteps
    for j = 1 : nModes
        if j == 1
            x_comb_IMM(i, :) = (mu(i, j) * x_IMM{i, j} )';
            P_comb_IMM(i, :, :) = mu(i, j) * (P_IMM{i, j} + (x_IMM{i, j} - (x_comb_IMM(i, :) )' ) * (x_IMM{i, j} - (x_comb_IMM(i, :) )' )' );
        else
            x_comb_IMM(i, :) = x_comb_IMM(i, :) + (mu(i, j) * x_IMM{i, j} )';
            P_comb_IMM(i, :, :) = shiftdim(P_comb_IMM(i, :, :) ) + mu(i, j) * (P_IMM{i, j} + (x_IMM{i, j} - (x_comb_IMM(i, :) )' ) * (x_IMM{i, j} - (x_comb_IMM(i, :) )' )' );
        end
    end
end

% step 6
GM = cell(nSteps, 1);
% construct Gaussian mixture
for i = 1 : nSteps
    mean = cell(nModes, 1);
    cov = cell(nModes, 1);
    for j = 1 : nModes
        mean{j} = x_IMM{i, j};
        cov{j} = P_IMM{i, j};
    end
    gm = func_constr_gm(mean, cov, (mu(i, :) )' );
    % log
    GM{i} = gm;
end

%% calculate the average estimated values of the ballistic coefficient
beta_est = zeros(nSteps, 1);
for i = 1 : nSteps
    for j = 1 : nModes
        if j == 1
            beta_est(i) = mu(i, j) * beta_IMM(j);
        else
            beta_est(i) = beta_est(i) + mu(i, j) * beta_IMM(j);
        end
    end
end

%% calculate MSE
error = x_tgt - x_comb_IMM;
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
    % average estimated values of the ballistic coefficient
    figure
    plot(T * (1 : nSteps), beta_est)
    title('average values of the estimated ballistic coefficient')
    xlabel('time (s)')
    ylabel('average values of the estimated ballistic coefficient (kg*m^-1*s^-2)')
    grid on
    % IMM mode probabilities
    figure
    if nModes == 2
        plot(T * (1 : nSteps), mu(:, 1), T * (1 : nSteps), mu(:, 2) )
    elseif nModes == 3
        plot(T * (1 : nSteps), mu(:, 1), T * (1 : nSteps), mu(:, 2), T * (1 : nSteps), mu(:, 3) )
    end
    title('IMM probabilities')
    xlabel('time (s)')
    ylabel('probabilities')
    legend('mode 1', 'mode 2')
    grid on
    % IMM
    figure
    plot(x_tgt(:, 1), x_tgt(:, 3) )
    hold on
    plot(x_comb_IMM(:, 1), x_comb_IMM(:, 3), '*')
    legend('true target trajectory', 'IMM')
    hold off
    xlabel('x (m)')
    ylabel('y (m)')
    grid on
    
end

cd('/Users/admin/Documents/GitHub/BTT_GCI_glint')
% log
logTime = sprintf('%s', datestr(now,30));
if ismac
    dataName = strcat('log/data', logTime);
else
    dataName = strcat('log\data', logTime);
end
save(dataName);