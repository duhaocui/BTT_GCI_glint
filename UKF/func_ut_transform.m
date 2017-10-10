function [m_tfr, P_tfr, P_cov_tfr] = func_ut_transform(m, P, func_handle, func_param, tr_param)
% input:
%           m - mean 
%           P - covariance
% func_handle - transformation function of the form g(x, param)
%  func_param - parameters of g, optional, default empty
%    tr_param - parameters of the unscented transform as:
%       kappa = tr_param{1}
% output:
%       m_tfr - transformed mean of y
%       P_trf - transformed covariance of y
%   P_cov_tfr - transformed cross-covariance of x and y 

% construct sigma points
[X, W] = func_constr_sigma(m, P, tr_param);

% calculate `Y` (propagate through the nonlinear function `g`)
if ischar(func_handle) || isa(func_handle,'function_handle')
    Y = [];
    for i = 1:size(X, 2)
        Y = [Y, feval(func_handle, X(:, i), func_param) ];
    end
else
    msg = 'ERROR occurred: ';
    error(msg)
end

n = length(m);
% calculate `m_tfr`
m_tfr = [];
for i = 1 : 2 * n + 1
    if i == 1
        m_tfr = W(i) * Y(:, i);
    else
        m_tfr = m_tfr + W(i) * Y(:, i);
    end
end
% calculate `P_tfr`
P_tfr = [];
for i = 1 : 2 * n + 1
    if i == 1
        P_tfr = W(i) * (Y(:, i) - m_tfr) * (Y(:, i) - m_tfr)';
    else
        P_tfr = P_tfr + W(i) * (Y(:, i) - m_tfr) * (Y(:, i) - m_tfr)';
    end
end
% calculate `P_cov_tfr`
P_cov_tfr = [];
for i = 1 : 2 * n + 1
    if i == 1
        P_cov_tfr = W(i) * (X(:, i) - m) * (Y(:, i) - m_tfr)';
    else
        P_cov_tfr = P_cov_tfr + W(i) * (X(:, i) - m) * (Y(:, i) - m_tfr)';
    end
end


