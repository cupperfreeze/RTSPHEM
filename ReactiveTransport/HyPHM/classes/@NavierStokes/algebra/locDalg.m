% def. of symbols
syms lam1 lam2 lam3 phi1 phi2 phi3 phi4 phi5 phi6 x y real
syms dxlam1 dxlam2 dxlam3 dxphi1 dxphi2 dxphi3 dxphi4 dxphi5 dxphi6 real
syms dylam1 dylam2 dylam3 dyphi1 dyphi2 dyphi3 dyphi4 dyphi5 dyphi6 real

% def. of vectors (for index in for-loops)
lam = [lam1, lam2, lam3];
dxlam = [dxlam1, dxlam2, dxlam3];
dylam = [dylam1, dylam2, dylam3];
phi = [phi1, phi2, phi3, phi4, phi5, phi6];
dxphi = [dxphi1, dxphi2, dxphi3, dxphi4, dxphi5, dxphi6];
dyphi = [dyphi1, dyphi2, dyphi3, dyphi4, dyphi5, dyphi6];

% def. of phi in terms of lambdas
for k = 1:3
    phi(k) = lam(k) * (2 * lam(k) - 1);
end
phi(4) = 4 * lam(2) * lam(3);
phi(5) = 4 * lam(1) * lam(3);
phi(6) = 4 * lam(1) * lam(2);

% def. of dxphi and dyphi in terms of lambdas
for k = 1:3
    dxphi(k) = dxlam(k) * (4 * lam(k) - 1);
    dyphi(k) = dylam(k) * (4 * lam(k) - 1);
end
dxphi(4) = 4 * (lam(2) * dxlam(3) + lam(3) * dxlam(2));
dyphi(4) = 4 * (lam(2) * dylam(3) + lam(3) * dylam(2));
dxphi(5) = 4 * (lam(1) * dxlam(3) + lam(3) * dxlam(1));
dyphi(5) = 4 * (lam(1) * dylam(3) + lam(3) * dylam(1));
dxphi(6) = 4 * (lam(1) * dxlam(2) + lam(2) * dxlam(1));
dyphi(6) = 4 * (lam(1) * dylam(2) + lam(2) * dylam(1));


Du = cell(6, 1);
Dv = cell(6, 1);

for j = 1:6


    Du{j} = sym(zeros(6));
    Dv{j} = sym(zeros(6));
    for k = 1:6
        for ell = 1:6
            Du{j}(ell, k) = phi(j) * dxphi(k) * phi(ell);
            Dv{j}(ell, k) = phi(j) * dyphi(k) * phi(ell);
        end
    end
    Du{j} = expand(Du{j}); % this is what can be integrated algebraically (on paper)
    Dv{j} = expand(Dv{j});


    % setting lambdas for REFERENCE TRIANGLE


    %   Du{j} = Du{j};   Dv{j} = Dv{j};

    Du{j} = subs(Du{j}, lam1, 1-x-y, 0);
    Du{j} = subs(Du{j}, lam2, x, 0);
    Du{j} = subs(Du{j}, lam3, y, 0);
    Du{j} = subs(Du{j}, dxlam1, -1, 0);
    Du{j} = subs(Du{j}, dxlam2, 1, 0);
    Du{j} = subs(Du{j}, dxlam3, 0, 0);
    Du{j} = subs(Du{j}, dylam1, -1, 0);
    Du{j} = subs(Du{j}, dylam2, 0, 0);
    Du{j} = subs(Du{j}, dylam3, 1, 0);

    Dv{j} = subs(Dv{j}, dylam1, -1, 0);
    Dv{j} = subs(Dv{j}, dylam2, 0, 0);
    Dv{j} = subs(Dv{j}, dylam3, 1, 0);
    Dv{j} = subs(Dv{j}, lam1, 1-x-y, 0);
    Dv{j} = subs(Dv{j}, lam2, x, 0);
    Dv{j} = subs(Dv{j}, lam3, y, 0);
    Dv{j} = subs(Dv{j}, dxlam1, -1, 0);
    Dv{j} = subs(Dv{j}, dxlam2, 1, 0);
    Dv{j} = subs(Dv{j}, dxlam3, 0, 0);

    Du{j} = int(Du{j}, y, 0, 1-x);
    Du{j} = int(Du{j}, x, 0, 1);

    Dv{j} = int(Dv{j}, y, 0, 1-x);
    Dv{j} = int(Dv{j}, x, 0, 1);


    display(Du{j})
    display(Dv{j})

end