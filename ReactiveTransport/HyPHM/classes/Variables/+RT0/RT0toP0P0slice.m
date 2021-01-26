%> @file RT0toP0P0slice.m Mapping @f$\mathbb{RT}_0(\mathcal{T}) \rightarrow \mathbb{P}_0(\mathcal{T})^2@f$.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> RT0toP0P0slice() evaluates the flux for each trangle of the grid
%> where the global @f$\mathbb{RT}_0(\mathcal{T})@f$ solution is given.  The cartesian flux on each
%> barycenter for each triangle is evaluated via getCartesianCoords().
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @param   g  Grid [ AbstractGrid ]
%> @param   fluxRT0  flux solution @f$\mathbb{RT}_0(\mathcal{T})@f$ basis  @f$[\#E\times1]@f$
%> @retval  flux   (flux) vector at each triangle  @f$[\#T\times2]@f$
%


function [fluxP0P0] = RT0toP0P0slice(g, fluxRT0)

assert(isa(g, 'AbstractGrid'))
RT0.checkdata(g, fluxRT0)

fluxP0P0 = zeros(g.numT, 2);

% for kT = 1 : g.numT
%   fluxP0P0(kT, :) = RT0.getCartesianCoords(g, fluxRT0(g.E0T(kT, :))', kT, g.baryT(kT,:)');
% end

coordV = permute(reshape(g.coordV0T(:, :, :), [g.numT, 3, 2]), [1, 3, 2]); % improved SG

for k = 1:3
    fluxP0P0 = fluxP0P0 + fluxRT0(g.E0T(:, k)) .* g.sigE0T(:, k) .* g.areaE(g.E0T(:, k)) .* (g.baryT(:, :) - coordV(:, :, k));
end
fluxP0P0 = 0.5 * fluxP0P0 ./ g.areaT;

return

end
