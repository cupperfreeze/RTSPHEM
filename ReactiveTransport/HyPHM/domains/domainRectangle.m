%> @file domainRectangle.m Generation of a meshed, rectangular bounded domain.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> [Grid] = domainRectangle(xmin, xmax, ymin, ymax, hmax)
%>
%> Rectangular bounded domain with lower left corner (xmin, ymin) and upper
%> right corner (xmax, ymax).
%>
%> @param  xmin   x-coordinate left border    [ scalar ]
%> @param  xmax   x-coordinate right border   [ scalar ]
%> @param  ymin   y-coordinate top border     [ scalar ]
%> @param  ymax   y-coordinate bottom border  [ scalar ]
%> @param  hmax   mesh fineness               [ scalar ]
%> @retval g      the grid in HyPHM format    [ Grid ]
%>
%> See also HyPHM, AbstractGrid, FoldedGrid, Grid
function g = domainRectangle(xmin, xmax, ymin, ymax, hmax)

assert(isscalar(xmin))
assert(isscalar(xmax))
assert(isscalar(ymin))
assert(isscalar(ymax))
assert(isscalar(hmax) && hmax > 0)

gd = [3; 4; xmin; xmax; xmax; xmin; ymin; ymin; ymax; ymax]; % geometry description
sf = 'rectangle'; % set formula
ns = names2ns('rectangle'); % name space

[p, e, t] = initmesh(decsg(gd, sf, ns), 'Hmax', hmax);

g = Grid(p, e, t);

end
