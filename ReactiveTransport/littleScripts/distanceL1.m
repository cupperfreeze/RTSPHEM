% auxiliary file
% Calculate distance by interpolation, routine adopted from HyPHM

%> @file +P0/distanceL2.m See Variable.distance.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @param grid Instance of Grid.
%> @param whdata Discrete P0-data @f$w_h@f$.
%> @param wfun Analytic function @f$w@f$, @c x -> w(x).

function ret = distanceL1(grid, whdata, wfun, rangeL, rangeR)

g = grid;
assert(isequal(size(whdata), [g.numT, 1]))

integrants = zeros(g.numT, 1);
for kT = 1:g.numT
    if g.baryT(kT, 1) > rangeL & g.baryT(kT, 1) < rangeR
        whfun = @(X) whdata(kT);
        intgrd = @(X) abs(wfun(X)-whfun(X)); %#ok<PFBNS> % (wh-w)^2
        integrants(kT) = intT(g, kT, intgrd, '34'); % int_T (wh-w)^2
    end
end
ret = sum(integrants); % sqrt sum_T int_T (wh-w)^2

end
