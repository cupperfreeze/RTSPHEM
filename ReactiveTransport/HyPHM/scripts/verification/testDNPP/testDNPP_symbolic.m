%> @file testDNPP_symbolic.m Prescription of an exact solution and computation of source terms and boundary coefficients via symbolic toolbox.

%% definition of symbolic variables
syms t x y real % variables
syms u v real % unknowns for x and y velocity
syms Pi real % symbolic pi

xvec = [x; y]; % vector of variables

%% abbreviations
cx = cos(Pi*x);
cy = cos(Pi*y);
sx = sin(Pi*x);
sy = sin(Pi*y);

%% prescibed exact solutions
u = -t * cx * sy;
v = t * sx * cy;
p = -1 / 4 * (cos(2 * Pi * x) + cos(2 * Pi * y));
c1 = t * cx;
c2 = t * sy;
phi = t / Pi^2 * (cx - sy);
% hence one obtains the flux unknowns
E = -grad(phi, xvec);
q1 = -grad(c1, xvec) + ([u; v] + E) * c1;
q2 = -grad(c2, xvec) + ([u; v] - E) * c2;

%% computation of the right-hand-sides/balancing terms of the equations % time derivatiev would be 'dt([u;v], t) '
DARCY = [u; v] + grad(p, xvec) - (c1 - c2) * E;

NP1 = dt(c1, t) + div(q1, xvec);
NP2 = dt(c2, t) + div(q2, xvec);

POISSON = div(E, xvec) - c1 + c2;

%% check the beautiful results
% pretty(DARCY)
DARCY_F = [t * (-cx * sy) + Pi / 2 * sin(2 * Pi * x) + t^2 / Pi * (sy - cx) * sx; ...
    t * (sx * cy) + Pi / 2 * sin(2 * Pi * y) + t^2 / Pi * (sy - cx) * cy];


NP1_F = (1 + Pi^2 * t) * cx + t^2 * (cx^2 - sx^2 - cx * sy + Pi * sx * cx * sy);
NP2_F = (1 + Pi^2 * t) * sy + t^2 * (sy^2 - cy^2 - cx * sy + Pi * sx * cy * cy);

printline(2, 'Validation of the evaluated right-hand sides')

pretty(expand(DARCY - DARCY_F))
expand(NP1-NP1_F)
expand(NP2-NP2_F)
expand(POISSON)

printline(1, 'Computation of the essential boundary conditions')
% u = subs(u, Pi, pi);    v = subs(v, Pi, pi);    p = subs(p, Pi, pi);
% c1 = subs(c1, Pi, pi); c2 = subs(c2, Pi, pi); phi = subs(phi, Pi, pi);
% q1 = subs(q1, Pi, pi);    q2 = subs(q2, Pi, pi);

% Darcy
printline(2, 'u_N south, y = 0')
pretty(dot(subs([u; v], y, 0), [0; -1]))
printline(2, 'u_N east, x = 1')
pretty(dot(subs([u; v], x, 1), [1; 0]))
printline(2, 'u_N north, y = 1')
pretty(dot(subs([u; v], y, 1), [0; 1]))
printline(2, 'u_N west, x = 0')
pretty(dot(subs([u; v], x, 0), [-1; 0]))

% Nernst-Planck c+
printline(2, 'c+_flux south, y = 0')
pretty(dot(subs(q1, y, 0), [0; -1]))
printline(2, 'c+_flux east, x = 1')
pretty(dot(subs(q1, x, 1), [1; 0]))
printline(2, 'c+_flux north, y = 1')
pretty(dot(subs(q1, y, 1), [0; 1]))
printline(2, 'c+_flux west, x = 0')
pretty(dot(subs(q1, x, 0), [-1; 0]))

% Nernst-Planck c-
printline(2, 'c-_flux south, y = 0')
pretty(dot(subs(q2, y, 0), [0; -1]))
printline(2, 'c-_flux east, x = 1')
pretty(dot(subs(q2, x, 1), [1; 0]))
printline(2, 'c-_flux north, y = 1')
pretty(dot(subs(q2, y, 1), [0; 1]))
printline(2, 'c-_flux west, x = 0')
pretty(dot(subs(q2, x, 0), [-1; 0]))

% Poisson
printline(2, 'E_N south, y = 0')
pretty(dot(subs(-grad(phi, xvec), y, 0), [0; -1]))
printline(2, 'E_N east, x = 1')
pretty(dot(subs(-grad(phi, xvec), x, 1), [1; 0]))
printline(2, 'E_N north, y = 1')
pretty(dot(subs(-grad(phi, xvec), y, 1), [0; 1]))
printline(2, 'E_N west, x = 0')
pretty(dot(subs(-grad(phi, xvec), x, 0), [-1; 0]))
