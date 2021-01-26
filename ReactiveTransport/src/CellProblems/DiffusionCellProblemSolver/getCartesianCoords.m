%> @file getCartesianCoords.m Determines the cartesian coordinates of a vector given in the local @f$\vec{RT}_0@f$-basis.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @param   g           Instance of grid [Grid].
%> @param   RT0Coords   coordinates of vector in @f$\vec{RT}_0@f$-basis [1x3].
%> @param   kT          number of triangle [scalar].
%> @param   xEval       point within the triangle where to evaluate the flux @f$\vec{q}@f$ [2x1].
%> @retval  q           (flux) vector at xEval in triangle [2x1].

function [q] = getCartesianCoords(g, RT0Coords, kT, xEval)

% assert(isa(g, 'AbstractGrid'), ...                      improved SG
%   'Wrong arguments in getCartesianCoords(): g.')
assert(isequal(size(RT0Coords), [1, 3]), ...
    'Wrong arguments in getCartesianCoords(): RT0Coords.')
assert(isequal(size(xEval), [2, 1]), ...
    'Wrong arguments in getCartesianCoords(): xEval.')

% COMPUTATION    q = sum alpha_i phi_i
coordV = reshape(g.coordV0T(kT, :, :), 3, 2)';
q = [0; 0];
for k = 1:3 % improved SG
    %coordV = reshape(g.coordV0T(kT, k, :),2,1);                             % old implementation using global vertex coordinates: coordV = g.coordV(g.V0T(kT, k), :)';
    q = q + RT0Coords(k) * g.sigE0T(kT, k) * g.areaE(g.E0T(kT, k)) * (xEval - coordV(:, k));
end
q = 0.5 * q / g.areaT(kT);

return

end
