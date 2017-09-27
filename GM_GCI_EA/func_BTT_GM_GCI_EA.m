function [mse1_px, mse1_py, beta1_est, mse2_px, mse2_py, beta2_est, mse_f_px, mse_f_py] = func_BTT_GM_GCI_EA(T, q, g, beta_tgt, F, G, nSensors, p_R, sigma_r, sigma_theta, nModes, beta1_IMM, beta2_IMM, DEMO_FLAG)

%% target
% parameters
nSteps = 100 / T;
Theta = [T^3 / 3, T^2 / 2; T^2 / 2, T];
Q = q * blkdiag(Theta, Theta);
% initial target state
x_0_tgt = [232 * 1e3; 2290 * cosd(190); 88 * 1e3; 2290 * sind(190)];
nStates = length(x_0_tgt);
x_tgt = zeros(nSteps, nStates);
aero_drag_force = zeros(nSteps, 2);
% update target state
for i = 1 : nSteps
    if i == 1
        x_prev = x_0_tgt;
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
Z = cell(2, 1);
Z{1} = zeros(nSteps, 2);
Z{2} = zeros(nSteps, 2);
% update measurement
for i = 1 : nSteps
    for j = 1 : nSensors
        x = x_tgt(i, 1);
        y = x_tgt(i, 3);
        v_r = normrnd(0, sigma_r(j) );
        v_theta = normrnd(0, sigma_theta(j) );
        x_R = p_R(j, 1);
        y_R = p_R(j, 2);
        z_r = sqrt( (x - x_R)^2 + (y - y_R)^2) + v_r;
        z_theta = atan2( (y - y_R), (x - x_R) ) + v_theta;
        % log
        Z{j}(i, 1) = z_r;
        Z{j}(i, 2) = z_theta;
    end
end

%% local filter: IMM; fusion rule: GM fusion
% parameters
R1 = diag([sigma_r(1)^2, sigma_theta(1)^2]);
R2 = diag([sigma_r(2)^2, sigma_theta(2)^2]);

p1 = zeros(nModes(1), nModes(1) );
if nModes(1) == 2
    p1(1, 1) = 0.9;
    p1(1, 2) = 1 - p1(1, 1);
    p1(2, 1) = p1(1, 2);
    p1(2, 2) = p1(1, 1);
elseif nModes(1) == 3
    p1 = [0.9, 0.05, 0.05;
        0.05, 0.9, 0.05;
        0.05, 0.05, 0.9];
end

p2 = zeros(nModes(2), nModes(2) );
if nModes(2) == 2
    p2(1, 1) = 0.9;
    p2(1, 2) = 1 - p2(1, 1);
    p2(2, 1) = p2(1, 2);
    p2(2, 2) = p2(1, 1);
elseif nModes(2) == 3
    p2 = [0.9, 0.05, 0.05;
        0.05, 0.9, 0.05;
        0.05, 0.05, 0.9];
end

% initial estimate
x_0_IMM = zeros(nStates, 1);
x_0_IMM(1) = 230 * 1e3;
x_0_IMM(2) = 2300 * cosd(190);
x_0_IMM(3) = 90 * 1e3;
x_0_IMM(4) = 2300 * sind(190);
P_0_IMM = diag([1000^2, 20^2, 1000^2, 20^2]);

mu1_0 = ones(nModes(1), 1) / nModes(1);
mu2_0 = ones(nModes(2), 1) / nModes(2);

% sensors
x1_IMM_prior = cell(nSteps, nModes(1) );
P1_IMM_prior = cell(nSteps, nModes(1) );
x1_IMM = cell(nSteps, nModes(1) );
P1_IMM = cell(nSteps, nModes(1) );
mu1 = zeros(nSteps, nModes(1) );

x2_IMM_prior = cell(nSteps, nModes(2) );
P2_IMM_prior = cell(nSteps, nModes(2) );
x2_IMM = cell(nSteps, nModes(2) );
P2_IMM = cell(nSteps, nModes(2) );
mu2 = zeros(nSteps, nModes(2) );

x1_comb_IMM = zeros(nSteps, nStates);
x2_comb_IMM = zeros(nSteps, nStates);
x_f_comb_IMM = zeros(nSteps, nStates);
P1_comb_IMM = zeros(nSteps, nStates, nStates);
P2_comb_IMM = zeros(nSteps, nStates, nStates);
P_f_comb_IMM = zeros(nSteps, nStates, nStates);
GM1 = cell(nSteps, 1);
GM2 = cell(nSteps, 1);
GM_f = cell(nSteps, 1);

% IMM
for i = 1 : nSteps
    % step 1
    if i == 1
        mu1_km1 = mu1_0;
        mu2_km1 = mu2_0;
    else
        mu1_km1 = (mu1(i - 1, :) )';
        mu2_km1 = (mu2(i - 1, :) )';
    end
    % calculate the normalization constants c1_bar, c2_bar
    c1_bar = zeros(nModes(1), 1);
    c2_bar = zeros(nModes(2), 1);
    
    for j = 1 : nModes(1)
        for k = 1 : nModes(1)
            if k == 1
                c1_bar(j) = p1(k, j) * mu1_km1(k);
            else
                c1_bar(j) = c1_bar(j) + p1(k, j) * mu1_km1(k);
            end
        end
    end
    
    for j = 1 : nModes(2)
        for k = 1 : nModes(2)
            if k == 1
                c2_bar(j) = p2(k, j) * mu2_km1(k);
            else
                c2_bar(j) = c2_bar(j) + p2(k, j) * mu2_km1(k);
            end
        end
    end
    
    
    % calculate the mixing probabilities
    mu1_cond_km1 = zeros(nModes(1), nModes(1) );
    mu2_cond_km1 = zeros(nModes(2), nModes(2) );
    
    for j = 1 : nModes(1)
        for k = 1 : nModes(1)
            mu1_cond_km1(j, k) = p1(j, k) * mu1_km1(j) / c1_bar(k);
        end
    end
    
    for j = 1 : nModes(2)
        for k = 1 : nModes(2)
            mu2_cond_km1(j, k) = p2(j, k) * mu2_km1(j) / c2_bar(k);
        end
    end
    
    % step 2
    % get previous estimate matched to mode j
    x1_km1 = cell(nModes(1), 1);
    P1_km1 = cell(nModes(1), 1);    
    x2_km1 = cell(nModes(2), 1);
    P2_km1 = cell(nModes(2), 1);
    
    for j = 1 : nModes(1)
        if i == 1
            x1_km1{j} = x_0_IMM;
            P1_km1{j} = P_0_IMM;
        else
            x1_km1{j} = x1_IMM{i - 1, j};
            P1_km1{j} = P1_IMM{i - 1, j};
        end
    end
    
    for j = 1 : nModes(2)
        if i == 1
            x2_km1{j} = x_0_IMM;
            P2_km1{j} = P_0_IMM;
        else
            x2_km1{j} = x2_IMM{i - 1, j};
            P2_km1{j} = P2_IMM{i - 1, j};
        end
    end
    
    
    % initial mixed initial condition
    x1_0_km1 = cell(nModes(1), 1);
    x2_0_km1 = cell(nModes(2), 1);
    
    for j = 1 : nModes(1)
        for k = 1 : nModes(1)
            if k == 1
                x1_0_km1{j} = mu1_cond_km1(k, j) * x1_km1{k};
            else
                x1_0_km1{j} = x1_0_km1{j} + mu1_cond_km1(k, j) * x1_km1{k};
            end
        end
    end
    
    for j = 1 : nModes(2)
        for k = 1 : nModes(2)
            if k == 1
                x2_0_km1{j} = mu2_cond_km1(k, j) * x2_km1{k};
            else
                x2_0_km1{j} = x2_0_km1{j} + mu2_cond_km1(k, j) * x2_km1{k};
            end
        end
    end
    
       
    P1_0_km1 = cell(nModes(1), 1);
    P2_0_km1 = cell(nModes(2), 1);
    
    for j = 1 : nModes(1)
        for k = 1 : nModes(1)
            if k == 1
                P1_0_km1{j} = mu1_cond_km1(k, j) * (P1_km1{k} + (x1_km1{k} - x1_0_km1{j} ) * (x1_km1{k} - x1_0_km1{j} )' );
            else
                P1_0_km1{j} = P1_0_km1{j} + mu1_cond_km1(k, j) * (P1_km1{k} + (x1_km1{k} - x1_0_km1{j} ) * (x1_km1{k} - x1_0_km1{j} )' );
            end
        end
    end
    
    for j = 1 : nModes(2)
        for k = 1 : nModes(2)
            if k == 1
                P2_0_km1{j} = mu2_cond_km1(k, j) * (P2_km1{k} + (x2_km1{k} - x2_0_km1{j} ) * (x2_km1{k} - x2_0_km1{j} )' );
            else
                P2_0_km1{j} = P2_0_km1{j} + mu2_cond_km1(k, j) * (P2_km1{k} + (x2_km1{k} - x2_0_km1{j} ) * (x2_km1{k} - x2_0_km1{j} )' );
            end
        end
    end
    
    
    % step 3
    % mode-matched filtering
    z1 = (Z{1}(i, :) )';
    z2 = (Z{2}(i, :) )';
    
    Lambda1 = zeros(nModes(1), 1);
    Lambda2 = zeros(nModes(2), 1);
    
    for j = 1 : nModes(1)
        % predict
        f1_J = func_calc_f_J(g, beta1_IMM(j), x1_0_km1{j});
        [x1_pred, P1_pred] = func_EKF_pred(x1_0_km1{j}, P1_0_km1{j}, F, G, g, f1_J, Q, beta1_IMM(j) );
        % update
        x1_R = p_R(1, 1);
        y1_R = p_R(1, 2);
        h1_J = func_calc_h_J(x1_R, y1_R, x1_pred);
        [x1_upd, P1_upd] = func_EKF_upd(x1_R, y1_R, x1_pred, P1_pred, z1, h1_J, R1);
        % calculate the likelihood functions matched to mode j
        z1_r = sqrt( (x1_pred(1) - x1_R)^2 + (x1_pred(3) - y1_R)^2);  
        z1_theta = atan2( (x1_pred(3) - y1_R), (x1_pred(1) - x1_R) );     
        mean1_prior = [z1_r; z1_theta];
        cov1_prior = h1_J * P1_pred * h1_J' + R1;
        Lambda1(j) = mvnpdf(z1', mean1_prior', cov1_prior);
        % log
        x1_IMM_prior{i, j} = x1_pred;
        P1_IMM_prior{i, j} = P1_pred;
        x1_IMM{i, j} = x1_upd;
        P1_IMM{i, j} = P1_upd;
    end
    
    for j = 1 : nModes(2)
        % predict
        f2_J = func_calc_f_J(g, beta2_IMM(j), x2_0_km1{j});
        [x2_pred, P2_pred] = func_EKF_pred(x2_0_km1{j}, P2_0_km1{j}, F, G, g, f2_J, Q, beta2_IMM(j) );
        % update
        x2_R = p_R(2, 1);
        y2_R = p_R(2, 2);
        h2_J = func_calc_h_J(x2_R, y2_R, x2_pred);
        [x2_upd, P2_upd] = func_EKF_upd(x2_R, y2_R, x2_pred, P2_pred, z2, h2_J, R2);
        % calculate the likelihood functions matched to mode j
        z2_r = sqrt( (x2_pred(1) - x2_R)^2 + (x2_pred(3) - y2_R)^2);       
        z2_theta = atan2( (x2_pred(3) - y2_R), (x2_pred(1) - x2_R) );      
        mean2_prior = [z2_r; z2_theta];
        cov2_prior = h2_J * P2_pred * h2_J' + R2;
        Lambda2(j) = mvnpdf(z2', mean2_prior', cov2_prior);
        % log
        x2_IMM_prior{i, j} = x2_pred;
        P2_IMM_prior{i, j} = P2_pred;     
        x2_IMM{i, j} = x2_upd;
        P2_IMM{i, j} = P2_upd;
    end
    
    
    
    % step 4
    % calculate the normalization constants c
    for j = 1 : nModes(1)
        if j == 1
            c1 = Lambda1(j) * c1_bar(j);
        else
            c1 = c1 + Lambda1(j) * c1_bar(j);
        end
    end
    
    for j = 1 : nModes(2)
        if j == 1
            c2 = Lambda2(j) * c2_bar(j);
        else
            c2 = c2 + Lambda2(j) * c2_bar(j);
        end
    end   
    
    % mode probability update
    for j = 1 : nModes(1)
        mu1(i, j) = Lambda1(j) * c1_bar(j) / c1;
    end
    
    for j = 1 : nModes(2)
        mu2(i, j) = Lambda2(j) * c2_bar(j) / c2;
    end
    
    % step 5
    % construct Gaussian mixture
    mean1 = cell(nModes(1), 1);
    cov1 = cell(nModes(1), 1);
    mean2 = cell(nModes(2), 1);
    cov2 = cell(nModes(2), 1);
    
    for j = 1 : nModes(1)
        mean1{j} = x1_IMM{i, j};
        cov1{j} = P1_IMM{i, j};
    end  
    
    for j = 1 : nModes(2)
        mean2{j} = x2_IMM{i, j};
        cov2{j} = P2_IMM{i, j};
    end
    
    
    gm1 = func_constr_gm(mean1, cov1, (mu1(i, :) )' );
    gm2 = func_constr_gm(mean2, cov2, (mu2(i, :) )' );
    % log
    GM1{i} = gm1;
    GM2{i} = gm2;
  
    % step 6
    % estimate and covariance combination
    [t1, t2] = func_GM2Gaussian(gm1);
    x1_comb_IMM(i, :) = t1';
    P1_comb_IMM(i, :, :) = t2;
    [t1, t2] = func_GM2Gaussian(gm2);
    x2_comb_IMM(i, :) = t1';
    P2_comb_IMM(i, :, :) = t2;

end

% save('data_GM.mat', 'GM1', 'GM2')

% fusion
delta = 0.1;
for i = 1 : nSteps
    gm1 = GM1{i};
    gm2 = GM2{i};
    % Gaussian mixture GCI
    % (1) Gaussian mixture fusion
    w_GCI = func_find_optimal_w(gm1, gm2, delta);
    gm_f = func_GM_GCI_EA(gm1, gm2, w_GCI);
    % (2) log
    GM_f{i} = gm_f;
    % estimate and covariance combination
    [t1, t2] = func_GM2Gaussian(gm_f);
    x_f_comb_IMM(i, :) = t1';
    P_f_comb_IMM(i, :, :) = t2;
    
end

%% calculate the average estimated values of the ballistic coefficient
beta1_est = zeros(nSteps, 1);
beta2_est = zeros(nSteps, 1);
for i = 1 : nSteps
    for j = 1 : nModes(1)
        if j == 1
            beta1_est(i) = mu1(i, j) * beta1_IMM(j);
        else
            beta1_est(i) = beta1_est(i) + mu1(i, j) * beta1_IMM(j);
        end
    end
end
for i = 1 : nSteps
    for j = 1 : nModes(2)
        if j == 1
            beta2_est(i) = mu2(i, j) * beta2_IMM(j);
        else
            beta2_est(i) = beta2_est(i) + mu2(i, j) * beta2_IMM(j);
        end
    end
end

%% calculate MSE
error1 = x_tgt - x1_comb_IMM;
error2 = x_tgt - x2_comb_IMM;
error_f = x_tgt - x_f_comb_IMM;
error1_px = error1(:, 1);
error2_px = error2(:, 1);
error_f_px = error_f(:, 1);
error1_py = error1(:, 3);
error2_py = error2(:, 3);
error_f_py = error_f(:, 3);

mse1_px = error1_px .* error1_px;
mse2_px = error2_px .* error2_px;
mse_f_px = error_f_px .* error_f_px;
mse1_py = error1_py .* error1_py;
mse2_py = error2_py .* error2_py;
mse_f_py = error_f_py .* error_f_py;


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
    subplot(211)
    plot(T * (1 : nSteps), aero_drag_force(:, 1) )
    xlabel('time (s)')
    ylabel('aerodynamic drag force in X (m/s^2)')
    grid on
    subplot(212)
    plot(T * (1 : nSteps), aero_drag_force(:, 2) )
    xlabel('time (s)')
    ylabel('aerodynamic drag force in Y (m/s^2)')
    grid on
    % range measurement
    figure
    subplot(211)
    plot(T * (1 : nSteps), Z{1}(:, 1) )
    title('range measurement')
    xlabel('time (s)')
    ylabel('range (m)')
    grid on
    % range measurement
    subplot(212)
    plot(T * (1 : nSteps), Z{2}(:, 1) )
    title('range measurement')
    xlabel('time (s)')
    ylabel('range (m)')
    grid on
    % bearing measurement
    figure
    subplot(221)
    plot(T * (1 : nSteps), Z{1}(:, 2))
    title('bearing measurement')
    xlabel('time (s)')
    ylabel('bearing (rad)')
    grid on
    subplot(222)
    plot(T * (1 : nSteps), rad2deg(Z{1}(:, 2) ) )
    title('bearing measurement')
    xlabel('time (s)')
    ylabel('bearing (deg)')
    grid on
    % bearing measurement
    subplot(223)
    plot(T * (1 : nSteps), Z{2}(:, 2))
    title('bearing measurement')
    xlabel('time (s)')
    ylabel('bearing (rad)')
    grid on
    subplot(224)
    plot(T * (1 : nSteps), rad2deg(Z{2}(:, 2) ) )
    title('bearing measurement')
    xlabel('time (s)')
    ylabel('bearing (deg)')
    grid on    
    % average estimated values of the ballistic coefficient
    figure
    subplot(211)
    plot(T * (1 : nSteps), beta1_est)
    title('average values of the estimated ballistic coefficient')
    xlabel('time (s)')
    ylabel('average values of the estimated ballistic coefficient (kg*m^-1*s^-2)')
    grid on
    subplot(212)
    plot(T * (1 : nSteps), beta2_est)
    title('average values of the estimated ballistic coefficient')
    xlabel('time (s)')
    ylabel('average values of the estimated ballistic coefficient (kg*m^-1*s^-2)')
    grid on
    % IMM mode probabilities
    figure
    subplot(211)
    if nModes(1) == 2
        plot(T * (1 : nSteps), mu1(:, 1), T * (1 : nSteps), mu1(:, 2) )
    elseif nModes(1) == 3
        plot(T * (1 : nSteps), mu1(:, 1), T * (1 : nSteps), mu1(:, 2), T * (1 : nSteps), mu1(:, 3) )
    end
    title('IMM probabilities')
    xlabel('time (s)')
    ylabel('probabilities')
    legend('mode 1', 'mode 2')
    grid on
    subplot(212)
    if nModes(2) == 2
        plot(T * (1 : nSteps), mu2(:, 1), T * (1 : nSteps), mu2(:, 2) )
    elseif nModes(2) == 3
        plot(T * (1 : nSteps), mu2(:, 1), T * (1 : nSteps), mu2(:, 2), T * (1 : nSteps), mu2(:, 3) )
    end
    title('IMM probabilities')
    xlabel('time (s)')
    ylabel('probabilities')
    legend('mode 1', 'mode 2')
    grid on
    % IMM
    figure
    subplot(211)
    plot(x_tgt(:, 1), x_tgt(:, 3) )
    hold on
    plot(x1_comb_IMM(:, 1), x1_comb_IMM(:, 3), '*')
    legend('true target trajectory', 'IMM')
    hold off
    xlabel('x (m)')
    ylabel('y (m)')
    grid on
    subplot(212)
    plot(x_tgt(:, 1), x_tgt(:, 3) )
    hold on
    plot(x2_comb_IMM(:, 1), x2_comb_IMM(:, 3), '*')
    legend('true target trajectory', 'IMM')
    hold off
    xlabel('x (m)')
    ylabel('y (m)')
    grid on
    
end