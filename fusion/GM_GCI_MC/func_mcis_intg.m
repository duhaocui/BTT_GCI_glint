function intg = func_mcis_intg(f1, f2, f_pi, w, N)

x = random(f_pi, N);

f_ij = zeros(N, 1);
for i = 1 : N
    t1 = pdf(f1, x(i, :) );
    t2 = pdf(f2, x(i, :) );
    f_ij(i) = (t1 / t2)^w;
end

for i = 1 : N
    t = f_ij(i) * pdf(f2, x(i, :) ) / pdf(f_pi, x(i, :) );
    if i == 1
        intg = t;
    else
        intg = intg + t;
    end
end
intg = intg / N