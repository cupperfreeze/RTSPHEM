%> @file +P0/distanceL2.m See Variable.distance.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @param grid Instance of Grid.
%> @param whdata Discrete P0-data @f$w_h@f$.
%> @param wfun Analytic function @f$w@f$, @c x -> w(x).

function ret = distanceL2(grid, whdata, wfun)

g = grid;
assert(isequal(size(whdata), [g.numT, 1]))

integrants = zeros(g.numT, 1);
parfor kT = 1:g.numT
    whfun = @(X) whdata(kT);
    intgrd = @(X) (wfun(X) - whfun(X))^2; %#ok<PFBNS> % (wh-w)^2
    integrants(kT) = intT(g, kT, intgrd, '612'); % int_T (wh-w)^2
end
ret = sqrt(sum(integrants)); % sqrt sum_T int_T (wh-w)^2

end
