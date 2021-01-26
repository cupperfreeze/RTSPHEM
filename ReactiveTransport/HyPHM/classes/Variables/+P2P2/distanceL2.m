%> @file +P2P2/distanceL2.m See Variable.distance.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> Depends on P2P2/localdata2fh.m.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @param grid Instance of Grid.
%> @param whdata Discrete P2P2-data @f$w_h@f$ @f$[\#V+\#E \times 2]@f$.
%> @param wfun Analytic function @f$\vec{w}@f$, @f$ \vec{x} \mapsto \vec{w}(\vec{x})@f$.

function ret = distanceL2(grid, whdata, wfun) %#ok<*PFBNS>

g = grid;
assert(isequal(size(whdata), [g.numV + g.numE, 2]))
integrants = zeros(g.numT, 1);

parfor kT = 1:g.numT
    % vector of function values on the 6 nodes of both components [6x2]
    coordsT = whdata([g.V0T(kT, :), g.numV + g.E0T(kT, :)], :);
    whfun = P2P2.localdata2fh(g, coordsT, kT);

    intgrd = @(X) (idx(wfun(X), 1) - idx(whfun(X), 1))^2 + ...
        (idx(wfun(X), 2) - idx(whfun(X), 2))^2; % ( ||wh-w||_2 )^2
    integrants(kT) = intT(g, kT, intgrd, '612'); % int_T (wh-w)^2
end

ret = sqrt(sum(integrants)); % sqrt sum_T int_T (wh-w)^2

end
