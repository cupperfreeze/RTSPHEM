%> @file alg.m Algebraic computation of local stiffness matrices.

function Href = alg

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


% %%%%%%%%%%
% %% locA %%
% %%%%%%%%%%
%
% A = sym(zeros(6));
% for k = 1 : 6
%   for ell = 1 : 6
%     A(k, ell) = phi(k) * phi(ell);
%   end
% end
% A = expand(A); % this is what can be integrated algebraically (on paper)
%
% % setting lambdas for REFERENCE TRIANGLE
% Aref = A;
% Aref = subs(Aref, lam1, 1-x-y); Aref = subs(Aref, lam2, x); Aref = subs(Aref, lam3, y);
% Aref = int(Aref, y, 0, 1-x);
% Aref = int(Aref, x, 0, 1);
%
%
% %%%%%%%%%%
% %% locB %%
% %%%%%%%%%%
%
% B1 = sym(zeros(6));
% for k = 1 : 6
%   for ell = 1 : 6
%     B1(ell, k) = dxphi(k) * dxphi(ell);
%   end
% end
% B1 = expand(B1); % this is what can be integrated algebraically (on paper)
%
% % setting lambdas for REFERENCE TRIANGLE
% B1ref = B1;
% B1ref = subs(B1ref, lam1, 1-x-y); B1ref = subs(B1ref, lam2, x); B1ref = subs(B1ref, lam3, y);
% B1ref = subs(B1ref, dxlam1, -1); B1ref = subs(B1ref, dxlam2, 1); B1ref = subs(B1ref, dxlam3, 0);
% B1ref = int(B1ref, y, 0, 1-x);
% B1ref = int(B1ref, x, 0, 1);
%
%
% B2 = sym(zeros(6));
% for k = 1 : 6
%   for ell = 1 : 6
%     B2(ell, k) = dyphi(k) * dyphi(ell);
%   end
% end
% B2 = expand(B2); % this is what can be integrated algebraically (on paper)
%
% % setting lambdas for REFERENCE TRIANGLE
% B2ref = B2;
% B2ref = subs(B2ref, lam1, 1-x-y); B2ref = subs(B2ref, lam2, x); B2ref = subs(B2ref, lam3, y);
% B2ref = subs(B2ref, dylam1, -1); B2ref = subs(B2ref, dylam2, 0); B2ref = subs(B2ref, dylam3, 1);
% B2ref = int(B2ref, y, 0, 1-x);
% B2ref = int(B2ref, x, 0, 1);
%
%
%
% B3 = sym(zeros(6));
% for k = 1 : 6
%   for ell = 1 : 6
%     B3(ell, k) = dyphi(k) * dxphi(ell) + dyphi(ell) * dxphi(k) ;
%   end
% end
% B3 = expand(B3); % this is what can be integrated algebraically (on paper)
%
% % setting lambdas for REFERENCE TRIANGLE
% B3ref = B3;
% B3ref = subs(B3ref, lam1, 1-x-y); B3ref = subs(B3ref, lam2, x); B3ref = subs(B3ref, lam3, y);
% B3ref = subs(B3ref, dxlam1, -1); B3ref = subs(B3ref, dxlam2, 1); B3ref = subs(B3ref, dxlam3, 0);
% B3ref = subs(B3ref, dylam1, -1); B3ref = subs(B3ref, dylam2, 0); B3ref = subs(B3ref, dylam3, 1);
% B3ref = int(B3ref, y, 0, 1-x);
% B3ref = int(B3ref, x, 0, 1);
%
%
%
% %%%%%%%%%%
% %% locC %%
% %%%%%%%%%%
%
% C1 = sym(zeros(6));
% for k = 1 : 6
%   for ell = 1 : 6
%     C1(k, ell) = dyphi(k) * dyphi(ell);
%   end
% end
% C1 = expand(C1); % this is what can be integrated algebraically (on paper)
%
% % setting lambdas for REFERENCE TRIANGLE
% C1ref = C1;
% C1ref = subs(C1ref, lam1, 1-x-y); C1ref = subs(C1ref, lam2, x); C1ref = subs(C1ref, lam3, y);
% C1ref = subs(C1ref, dylam1, -1); C1ref = subs(C1ref, dylam2, 0); C1ref = subs(C1ref, dylam3, 1);
% C1ref = int(C1ref, y, 0, 1-x);
% C1ref = int(C1ref, x, 0, 1);
%
% C2 = sym(zeros(6));
% for k = 1 : 6
%   for ell = 1 : 6
%     C2(k, ell) = dxphi(k) * dxphi(ell);
%   end
% end
% C2 = expand(C2); % this is what can be integrated algebraically (on paper)
%
% % setting lambdas for REFERENCE TRIANGLE
% C2ref = C2;
% C2ref = subs(C2ref, lam1, 1-x-y); C2ref = subs(C2ref, lam2, x); C2ref = subs(C2ref, lam3, y);
% C2ref = subs(C2ref, dxlam1, -1); C2ref = subs(C2ref, dxlam2, 1); C2ref = subs(C2ref, dxlam3, 0);
% C2ref = int(C2ref, y, 0, 1-x);
% C2ref = int(C2ref, x, 0, 1);
%
% C3 = sym(zeros(6));
% for k = 1 : 6
%   for ell = 1 : 6
%     C3(k, ell) = dxphi(k) * dyphi(ell) + dxphi(ell) * dyphi(k);
%   end
% end
% C3 = expand(C3); % this is what can be integrated algebraically (on paper)
%
% % setting lambdas for REFERENCE TRIANGLE
% C3ref = C3;
% C3ref = subs(C3ref, lam1, 1-x-y); C3ref = subs(C3ref, lam2, x); C3ref = subs(C3ref, lam3, y);
% C3ref = subs(C3ref, dxlam1, -1); C3ref = subs(C3ref, dxlam2, 1); C3ref = subs(C3ref, dxlam3, 0);
% C3ref = subs(C3ref, dylam1, -1); C3ref = subs(C3ref, dylam2, 0); C3ref = subs(C3ref, dylam3, 1);
% C3ref = int(C3ref, y, 0, 1-x);
% C3ref = int(C3ref, x, 0, 1);
%
%
% %%%%%%%%%%
% %% locF %%
% %%%%%%%%%%
%
% F1 = sym(zeros(3, 6));
% for k = 1 : 3
%   for ell = 1 : 6
%     F1(k, ell) = lam(k) * dxphi(ell);
%   end
% end
% F1 = expand(F1); % this is what can be integrated algebraically (on paper)
%
% % setting lambdas for REFERENCE TRIANGLE
% F1 = subs(F1, lam1, 1-x-y); F1 = subs(F1, lam2, x); F1 = subs(F1, lam3, y);
% F1 = subs(F1, dxlam1, -1); F1 = subs(F1, dxlam2, 1); F1 = subs(F1, dxlam3, 0);
% F1 = int(F1, y, 0, 1-x);
% F1 = int(F1, x, 0, 1);
%
% F2 = sym(zeros(3, 6));
% for k = 1 : 3
%   for ell = 1 : 6
%     F2(k, ell) = lam(k) * dyphi(ell);
%   end
% end
% F2 = expand(F2); % this is what can be integrated algebraically (on paper)
%
% % setting lambdas for REFERENCE TRIANGLE
% F2 = subs(F2, lam1, 1-x-y); F2 = subs(F2, lam2, x); F2 = subs(F2, lam3, y);
% F2 = subs(F2, dylam1, -1); F2 = subs(F2, dylam2, 0); F2 = subs(F2, dylam3, 1);
% F2 = int(F2, y, 0, 1-x);
% F2 = int(F2, x, 0, 1);


%%%%%%%%%%

%% locH %%
%%%%%%%%%%

% Raviart Thomas basis
syms s1 s2 s3 psi1x psi1y psi2x psi2y psi3x psi3y J11 J12 J21 J22 di11 di12 di21 di22 real
syms arE1 arE2 arE3 arEr1 arEr2 arEr3 real % areas of edges and reference edges
% def. of vectors (for index in for-loops)
psi = [psi1x, psi2x, psi3x; ...
    psi1y, psi2y, psi3y];
sig = [s1, s2, s3];
arE = [arE1, arE2, arE3];
arEr = [arEr1, arEr2, arEr3];
J = [J11, J12; J21, J22]; % Jacobian
DI = [di11, di12; di21, di22]; % inverse diffusion tensor

H = sym(zeros(3));
for k = 1:3
    for ell = 1:3
        H(k, ell) = (sig(k) * (arE(k) / arEr(k)) * (DI * J) * psi(:, k))' * (sig(ell) * (arE(ell) / arEr(ell)) * J * psi(:, ell));
    end
end
H = expand(H); % this is what can be integrated algebraically (on paper)

% setting values for REFERENCE TRIANGLE
Href = H;


Href = subs(Href, psi1x, sqrt(2)*x);
Href = subs(Href, psi1y, sqrt(2)*y);
Href = subs(Href, psi2x, x-1);
Href = subs(Href, psi2y, y);
Href = subs(Href, psi3x, x);
Href = subs(Href, psi3y, y-1);


Href = int(Href, y, 0, 1-x);
Href = int(Href, x, 0, 1);


end % alg

% EXACT INTEGRATION
% returns int_T lam1^a*lam2^b*lam3^c / 2 / |T|
function ret = intex(a, b, c)
str = sprintf('%d!*%d!*%d!/(%d+%d+%d+2)!', a, b, c, a, b, c);
ret = sym(str);
end
