%> @file +P0/distanceL1.m See Variable.distance.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @param grid Instance of Grid.
%> @param whdata Discrete P0-data @f$w_h@f$.
%> @param wfun Analytic function @f$w@f$, @c x -> w(x).

function ret = distanceL1(grid, whdata, wfun)

g = grid;
assert(isequal(size(whdata), [g.numT, 1]))

integrants = zeros(g.numT, 1);
parfor kT = 1:g.numT
    whfun = @(X) whdata(kT);
    intgrd = @(X) abs(wfun(X)-whfun(X));
    integrants(kT) = intT(g, kT, intgrd, '612');
end
ret = sum(integrants);

end