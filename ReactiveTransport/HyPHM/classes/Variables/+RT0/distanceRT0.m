%> @file +RT0/distanceRT0.m See Variable.distance.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @param grid Instance of Grid.
%> @param whdata Discrete @f$\vec{RT}_0@f$-data for @f$\vec{w}_h@f$ @f$[\#E\times 1]@f$.
%> @param wfun Analytic function @f$\vec{w}@f$, @f$ \vec{x} \mapsto \vec{w}(\vec{x})@f$.

function ret = distanceRT0(grid, whdata, wfun) %#ok<*PFBNS>

g = grid;
assert(isequal(size(whdata), [g.numE, 1]))

integrants = zeros(g.numT, 1);
for kT = 1:g.numT
    for kE = g.E0T(kT, :)
        baryE = g.baryE(kE, :)';
        nuE = g.nuE(kE, :)';
        integrants(kT) = integrants(kT) + g.areaT(kT) / 3 * (abs(dot(wfun(baryE), nuE) - whdata(kE)))^2;
    end
end
ret = sqrt(sum(integrants)); % sqrt sum_T int_T (wh-w)^2

end
