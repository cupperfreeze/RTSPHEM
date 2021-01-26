% tol:	1E-6
% tau:  1/4


h = [1 / 2, 1 / 4, 1 / 8, 1 / 16, 1 / 32, 1 / 64, 1 / 128, 1 / 256];

testDNPP_printresults % initialize Err

u = Err(1, :);
p = Err(2, :);
q1 = Err(3, :);
c1 = Err(4, :);
q2 = Err(5, :);
c2 = Err(6, :);
E = Err(7, :);
phi = Err(8, :);


loglog(h, u, 'ko-')
hold on
loglog(h, p, 'ks-')
loglog(h, q1, 'k>-')
loglog(h, c1, 'kv-')
loglog(h, q2, 'k<-')
loglog(h, c2, 'k^-')
loglog(h, E, 'kd-')
loglog(h, phi, 'kh-')


xlim([0.0001, 1])

legend('||u_h - u||~~~~~~~~~~~~~~~~~~~~~~~~', ...
    '||p_h - p||', ...
    '||j^{+}_h - j^{+}||', ...
    '||c^{+}_h - c^{+}||', ...
    '||j^{-}_h - j^{-}||', ...
    '||c^{-}_h - c^{-}||', ...
    '||E_h - E||', ...
    '||\phi_h - \phi||')

c = 7;

% triangle with slope 1
line([h(end -5), h(end -6)], c*[c1(end -1), c1(end -2)]/c1(end -1), 'color', [0, 0, 0])
line([h(end -5), h(end -6)], c*[c1(end -1), c1(end -1)]/c1(end -1), 'color', [0, 0, 0])
line([h(end -6), h(end -6)], c*[c1(end -1), c1(end -2)]/c1(end -1), 'color', [0, 0, 0])
text(exp((log(h(end -5))+log(h(end -6))) / 2), (c - 4), '1')

matlabfrag('MMS_testDNPP')
