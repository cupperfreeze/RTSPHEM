% locbUsource performs the partial assembly of a local assembly vector
% under consideration of the global orientation.  This function may be
% called with Vold instead of Uold to assemble bVsource.
%
% Input
%   d        [ NavierStokes ]  data required by solver
%   kT       [ scalar     ]  number of current element
%   compNo   {1, 2}          component number, see below
%
% Output
%   stiff    [6 x 1]             local assembly vector <Fx, phi_N'>_kT
%     or <Fy, phi_N'>_kT for compNp equal to 1 or 2, respectively.
%     The function phi returns all 6 phi_N.
%
%   stiff = <Fx, phi_N'>_kT,  N' = 1,...,6    or
%   stiff = <Fy, phi_N'>_kT,  N' = 1,...,6
%
%   We choose 7 quadrature points to compute the integral as follows
%   (cf. Zhangxin Chen 2005, p. 50):
%
%     int f = |T|* (9/20*f(center) + 1/20*sum f(verts) + 2/15*sum f(baryE))
%
%
%                       Q3
%                                Q4
%                            Q7          Q2
%                      Q5
%                              Q6
%
%                     Q1
%
% where the quadrature point Q7 is approximated by taking the mean of the
% vertex values.
%
% See also phi, NavierStokes.computeLevel
%
% Copyright 2009, 2010 F. Frank, Chair of Applied Mathematics I,
%           Department of Mathematics, University of Erlangen-Nuremberg,
%           G E R M A N Y

function stiff = locbUsource(d, kT, time, compNo)

assert(isa(d, 'NavierStokes'))
assert(ismember(compNo, 1:2))
g = d.grid;

% INITIALIZATION
stiffV = zeros(6, 1); % the sum over vertices (Q1-Q3)
stiffE = zeros(6, 1); % the sum over edges    (Q4-Q6)
% stiffC :            % the center term       (Q7)

% QUADRATURE POINTS (3x VERTEX)
quadPoints = squeeze(g.coordV0T(kT, :, :))';
for k = 1:3
    phii = phi(g, kT, quadPoints(:, k));
    fvec = d.evalF(time, kT, k);
    fcomp = fvec(compNo);
    stiffV = stiffV + fcomp * phii;
end

% QUADRATURE POINTS (3x EDGE MIDPOINT)
%         | Q1(1)  ...  Q3(1) |
%         | Q1(2)  ...  Q3(2) |
quadPoints = squeeze(g.baryE0T(kT, :, :))';
% Fetch barycenters _related to triangle_.  This is important if grid is
% PeriodicGrid.
for k = 1:3
    phii = phi(g, kT, quadPoints(:, k));
    fvec = d.evalF(time, kT, k+3);
    fcomp = fvec(compNo);
    stiffE = stiffE + fcomp * phii;
end

% QUADRATURE POINT (1x CENTER)
phii = phi(g, kT, g.baryT(kT, :)');
fvec1 = d.evalF(time, kT, 1);
fcomp1 = fvec1(compNo);
fvec2 = d.evalF(time, kT, 2);
fcomp2 = fvec2(compNo);
fvec3 = d.evalF(time, kT, 3);
fcomp3 = fvec3(compNo);
fcomp = 1 / 3 * (fcomp1 + fcomp2 + fcomp3);
stiffC = fcomp * phii;

% int f = |T|* (9/20*f(center) + 1/20*sum f(verts) + 2/15*sum f(baryE))
stiff = g.areaT(kT) * (9 / 20 * stiffC + 1 / 20 * stiffV + 2 / 15 * stiffE);

end
