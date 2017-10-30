function X_out = func_ut_transform(X_in, func_handle, func_param)

% unscented transform
if ischar(func_handle) || isa(func_handle,'function_handle')
    X_out = [];
    for i = 1:size(X_in, 2)
        X_out = [X_out, feval(func_handle, X_in(:, i), func_param) ];
    end
else
    msg = 'ERROR occurred: ';
    error(msg)
end

% calculate `P_cov_tfr`
P_cov_out = [];
for i = 1 : L
    if i == 1
        P_cov_out = W_c(i) * (X_out(:, i) - m_out) * (X_out(:, i) - m_tfr)';
    else
        P_cov_out = P_cov_out + W_c(i) * (X(:, i) - m_out) * (X_out(:, i) - m_tfr)';
    end
end


