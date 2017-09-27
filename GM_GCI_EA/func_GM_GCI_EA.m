function gm_f = func_GM_GCI_EA(gm1, gm2, w)

% read data
w1 = gm1.w;
w2 = gm2.w;

x1 = gm1.x;
x2 = gm2.x;

P1 = gm1.P;
P2 = gm2.P;

N1 = length(w1);
N2 = length(w2);

% calculate P_f
P_f = cell(N1, N2);
for i = 1 : N1
    for j = 1 : N2
        inv_P_f = w * inv(P1{i} ) + (1 - w) * inv(P2{j} );
        P_f{i, j} = inv(inv_P_f);
    end
end

% calculate x_f
x_f = cell(N1, N2);
for i = 1 : N1
    for j = 1 : N2
        x_f{i, j} = P_f{i, j} * (w * inv(P1{i} ) * x1{i} + (1 - w) * inv(P2{j} ) * x2{j} );
    end
end

% calculate d
assert(length(x1{1} ) == length(x2{1} ) )
d = length(x1{1} );

% calcualte zeta1
zeta1 = zeros(N1, 1);
for i = 1 : N1
    zeta1(i) = - (d * log(2 * pi) - log(det(w * inv(P1{i} ) ) ) + w * (x1{i})' * (inv(P1{i} ) )' * x1{i} ) / 2;
end

% calcualte zeta2
zeta2 = zeros(N2, 1);
for i = 1 : N2
    zeta2(i) = - (d * log(2 * pi) - log(det( (1 - w) * inv(P2{i} ) ) ) + (1 - w) * (x2{i})' * (inv(P2{i} ) )' * x2{i} ) / 2;
end

% calculate zeta_f
zeta_f = zeros(N1, N2);
for i = 1 : N1
    for j = 1 : N2
        zeta_f(i, j) = - (d * log(2 * pi) - log(det(inv(P_f{i, j} ) ) ) + (x_f{i, j})' * (inv(P_f{i, j} ) )' * x_f{i, j} ) / 2;
    end
end

% calculate c_f
c_f = zeros(N1, N2);
for i = 1 : N1
    for j = 1 : N2
        t = zeta1(i) + zeta2(j) - zeta_f(i, j);
        c_f(i, j) = exp(t);
    end
end

% calculate w_f
w_f = zeros(N1, N2);
for i = 1 : N1
    for j = 1 : N2
        w_f(i, j) = (w1(i) )^(w) * (w2(j) )^(1 - w) * c_f(i, j);
    end
end

% remove the imaginary part of w_f
% w_f = real(w_f);

% sort w_f, x_f and P_f 
w_f_sorted = zeros(N1 * N2, 1);
x_f_sorted = cell(N1 * N2, 1);
P_f_sorted = cell(N1 * N2, 1);
for i = 1 : N1
    for j = 1 : N2
        t = j + (i - 1) * N2;
        w_f_sorted(t) = w_f(i, j);
        x_f_sorted{t} = x_f{i, j};
        P_f_sorted{t} = P_f{i, j};
    end
end

w_f_sorted = w_f_sorted / sum(w_f_sorted);

% for i = 1 : length(w_f_sorted)
%     if ( w_f_sorted(i) < 1e-4)
%         w_f_sorted(i) = [];
%         x_f_sorted(i) = [];
%         P_f_sorted(i) = [];
%     end
% end

% construct GM
gm_f = func_constr_gm(x_f_sorted, P_f_sorted, w_f_sorted);
         



