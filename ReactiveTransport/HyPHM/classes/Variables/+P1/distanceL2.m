%> @file +P1/distanceL2.m See Variable.distance.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> Depends on P1.localdata2fh.m.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @param grid Instance of Grid.
%> @param whdata Discrete P1-data @f$w_h@f$.
%> @param wfun Analytic function @f$w@f$, @c x -> w(x).

function ret = distanceL2(grid, whdata, wfun)

g = grid;

integrants = zeros(g.numT, 1);
parfor kT = 1:g.numT
    coordsT = whdata(g.V0T(kT, :)); % vector of function values on the 3 nodes
    whfun = P1.localdata2fh(g, coordsT, kT);
    intgrd = @(X) (wfun(X) - whfun(X))^2; % (wh-w)^2
    integrants(kT) = intT(g, kT, intgrd, '612'); % int_T (wh-w)^2
end
ret = sqrt(sum(integrants)); % sqrt sum_T int_T (wh-w)^2
end
