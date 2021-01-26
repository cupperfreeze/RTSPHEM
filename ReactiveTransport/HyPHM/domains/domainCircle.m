%> @file domainCircle.m Generation of a meshed domain with circular boundary.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> [Grid] = domainCircle(xcen, ycen, rad, hmax)
%>
%> Generation of a meshed domain with circular boundary around the center
%> (@c xcen, @c ycen) with radius @c rad and mesh size @hmax.
%>
%> @param  xcen   x coordinate of the circle's center  [ scalar ]
%> @param  ycen   y coordinate of the circle's center  [ scalar ]
%> @param  rad    radius of the circle                 [ scalar ]
%> @param  hmax   mesh fineness                        [ scalar ]
%> @retval g          [ Grid ]    the grid in HyPHM format
%>
%> See also HyPHM, AbstractGrid, FoldedGrid, Grid

function g = domainCircle(xcen, ycen, rad, hmax)

assert(isscalar(xcen))
assert(isscalar(ycen))
assert(isscalar(rad) && rad > 0)
assert(isscalar(hmax) && rad > 0)

gd = [1; xcen; ycen; rad]; % geometry description
[p, e, t] = initmesh(decsg(gd), 'Hmax', hmax); % 'doc decsg' for details

g = Grid(p, e, t);

end
