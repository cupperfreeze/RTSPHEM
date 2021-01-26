%> @file +RT0/distanceL2.m See Variable.distance.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> Depends on RT0/localdata2fh.m.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @param grid Instance of Grid.
%> @param whdata Discrete @f$\vec{RT}_0@f$-data for @f$\vec{w}_h@f$ @f$[\#E\times 1]@f$.
%> @param wfun Analytic function @f$\vec{w}@f$, @f$ \vec{x} \mapsto \vec{w}(\vec{x})@f$.

function ret = distanceL2(grid, whdata, wfun) %#ok<*PFBNS>

g = grid;
assert(isequal(size(whdata), [g.numE, 1]))
integrants = zeros(g.numT, 1);

parfor kT = 1:g.numT
    % vector of function values on the 3 nodes of both components [3x1]
    coordsT = whdata(g.E0T(kT, :));
    whfun = RT0.localdata2fh(g, coordsT, kT); % T -> IR^2

    intgrd = @(X) (idx(wfun(X), 1) - idx(whfun(X), 1))^2 + ...
        (idx(wfun(X), 2) - idx(whfun(X), 2))^2; % ( ||wh-w||_2 )^2
    integrants(kT) = intT(g, kT, intgrd, '612'); % int_T (wh-w)^2
end

ret = sqrt(sum(integrants)); % sqrt sum_T int_T (wh-w)^2

end
