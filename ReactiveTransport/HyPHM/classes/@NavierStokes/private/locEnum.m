% locE performs the assembly of a local assembly matrix
% under consideration of the global orientation.
%
% Input
%   d        [ NavierStokes ]  data required by solver
%   kT       [ 1   ]             number of current element
%
% Output
%   stiff    [6 x 6]             local assembly matrix
%     The function phi returns all 6 phi_N, i.e. we have to build the
%     tensor product to obtain the local 6x6 matrix.
%
%   stiff(N',N) = <dxphi_N dyphi_N'>,  N, N' = 1,...,6
%
%   The integrant is quadratic, thus, we choose 3 quadrature points to
%   compute the integral exactly as follows (cf. Zhangxin Chen 2005, p. 50):
%
%     int f = |T|/3  *  sum f(baryE)
%
% See also phi, NavierStokes
%
% Copyright 2009, 2010 F. Frank, Chair of Applied Mathematics I,
%           Department of Mathematics, University of Erlangen-Nuremberg,
%           G E R M A N Y

function stiff = locEnum(d, kT)

g = d.grid;

% INITIALIZATION
stiff = zeros(6);

% QUADRATURE POINTS (3x EDGE MIDPOINT)
%         | Q1(1)  ...  Q3(1) |
%         | Q1(2)  ...  Q3(2) |
quadPoints = squeeze(g.baryE0T(kT, :, :))';
% Fetch barycenters _related to triangle_.  This is important if grid is
% PeriodicGrid.

% LOOP OVER QUADRATURE POINTS
for kQ = 1:3
    xQ = quadPoints(:, kQ);
    dxphi = phi(g, kT, xQ, 'x');
    dyphi = phi(g, kT, xQ, 'y');

    stiff = stiff + dxphi * dyphi';
end

stiff = g.areaT(kT) / 3 * stiff';

end

% Remark: Index switch since Aloc(j,i) = int(f(i) * f(j))