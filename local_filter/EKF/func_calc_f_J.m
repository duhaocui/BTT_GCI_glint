function f_J = func_calc_f_J(g, beta, X)

x = X(1);
dx = X(2);
y = X(3);
dy = X(4);

[rho, c1, c2] = func_calc_rho(y);

f_J = zeros(2, length(X) );

f_J(1, 1) = 0;
f_J(1, 2) = -g * rho / (2 * beta) * (2 * dx^2 + dy^2) / sqrt(dx^2 + dy^2);
f_J(1, 3) = g * c2 * rho / (2 * beta) * dx * sqrt(dx^2 + dy^2);
f_J(1, 4) = -g * rho / (2 * beta) * dx * dy / sqrt(dx^2 + dy^2);

f_J(2, 1) = 0;
f_J(2, 2) = -g * rho / (2 * beta) * dx * dy / sqrt(dx^2 + dy^2);
f_J(2, 3) = g * c2 * rho / (2 * beta) * dy * sqrt(dx^2 + dy^2);
f_J(2, 4) = -g * rho / (2 * beta) * (dx + 2 * dy) / sqrt(dx^2 + dy^2);
 