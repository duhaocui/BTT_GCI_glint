function intg = func_mc_intg(f1, f2, w, N)

x = random(f2, N);

f_ij = zeros(N, 1);
for i = 1 : N
    t1 = pdf(f1, x(i, :) );
    t2 = pdf(f2, x(i, :) );
    f_ij(i) = (t1 / t2)^w;
end

for i = 1 : N
    t = f_ij(i);
    if i == 1
        intg = t;
    else
        intg = intg + t;
    end
end
intg = intg / N;