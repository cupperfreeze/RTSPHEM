syms lam1 lam2 lam3 phi1 phi2 phi3 phi4 phi5 phi6

lam = [lam1, lam2, lam3];
phi = [phi1, phi2, phi3, phi4, phi5, phi6];

for k = 1:3
    phi(k) = lam(k) * (2 * lam(k) - 1);
end
phi(4) = 4 * lam(2) * lam(3);
phi(5) = 4 * lam(1) * lam(3);
phi(6) = 4 * lam(1) * lam(2);

A = sym(zeros(6));

for k = 1:6
    for ell = 1:6
        A(k, ell) = phi(k) * phi(ell);
    end
end

A = expand(A);
