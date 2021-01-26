%> @file +P1/intmean.m See Variable.mean.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%> @param grid Instance of Grid.
%> @param whdata Discrete P1-data @f$w_h@f$.

function ret = intmean(grid, whdata)

g = grid;

integrals = zeros(g.numT, 1);
parfor kT = 1:g.numT
    integrals(kT) = g.areaT(kT) / 3 * sum(whdata(g.V0T(kT, :)));
end
ret = sum(integrals) / sum(g.areaT);
end
