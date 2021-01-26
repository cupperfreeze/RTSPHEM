%> @file +P0/distanceDOF.m See Variable.distance.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @param grid Instance of Grid.
%> @param whdata Discrete P0-data @f$w_h@f$.
%> @param wfun Analytic function @f$w@f$, @c x -> w(x).

function ret = distanceP0(grid, whdata, wfun)

g = grid;
assert(isequal(size(whdata), [g.numT, 1]))

integrants = zeros(g.numT, 1);
parfor kT = 1:g.numT
    integrants(kT) = g.areaT(kT) * (abs(wfun(g.baryT(kT, :)') - whdata(kT)))^2; %#ok<*PFBNS>
end
ret = sqrt(sum(integrants)); % sqrt sum_T int_T (wh-w)^2

end
