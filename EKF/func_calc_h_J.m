function h_J = func_calc_h_J(x_R, y_R, X)

x = X(1);
dx = X(2);
y = X(3);
dy = X(4);

h_J = zeros(2, length(X) );

h_J(1, 1) = (x - x_R) / sqrt( (x - x_R)^2 + (y - y_R)^2);
h_J(1, 2) = 0;
h_J(1, 3) = (y - y_R) / sqrt( (x - x_R)^2 + (y - y_R)^2);
h_J(1, 4) = 0;

h_J(2, 1) = - (y - y_R) / ( (x - x_R)^2 + (y - y_R)^2);
h_J(2, 2) = 0;
h_J(2, 3) = (x - x_R) / ( (x - x_R)^2 + (y - y_R)^2);
h_J(2, 4) = 0;


