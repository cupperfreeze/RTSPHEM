%> @file domainCellY.m Generation of a unit square-bounded domain with inner elliptic hole.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @brief Generation of a unit square-bounded domain with inner elliptic hole.
%>
%> @code [Grid] = domainCellY(a, b, phi, hmax) @endcode
%>
%> Unit square-bounded domain @f$[0,1]\times[0,1]@f$ with inner elliptic obstacle/hole.  The ellipse has
%> axes of length a and b and is rotated by phi.  The grid is not yet
%> folded.
%> @n @n
%> If You want to use periodic boundaries, do
%> @code
%>   fg = FoldedGrid(domainCellY(a, b, phi, hmax))
%> @endcode
%>
%> See also HyPHM, AbstractGrid, FoldedGrid, Grid
%>
%> @param    a     Length of axis a [scalar].
%> @param    b     Length of axis b [scalar].
%> @param    phi   Rotation angle of obstacle [scalar].
%> @param    hmax  Fineness of the grid, required by mesh-generator [scalar].
%> @param    cx    optional: x-coordinate of the center of Y [scalar] (<i>default:</i> cx = 0.5).
%> @param    cy    ... y-coordinate ... (<i>default:</i> cy = 0.5)


%> @retval   g     Instance of the Grid.
function g = domainCellY(a, b, phi, hmax, cx, cy)

assert(isscalar(a))
assert(isscalar(b))
assert(isscalar(phi))
assert(a < 1 && a > 0)
assert(b < 1 && b > 0)

if ~exist('cx', 'var')
    cx = 0.5;
    cy = 0.5;
else
    assert(isscalar(cx))
    assert(isscalar(cy))
end


gd = [3.0000, 4.0000; ... % geometry description
    4.0000, cx; ...
    cx - 0.5, cy; ...
    cx + 0.5, a / 2; ...
    cx + 0.5, b / 2; ...
    cx - 0.5, phi; ...
    cy - 0.5, 0; ...
    cy - 0.5, 0; ...
    cy + 0.5, 0; ...
    cy + 0.5, 0];

sf = 'square-ellipse'; % set formula
ns = names2ns('square', 'ellipse'); % name space

[p, e, t] = initmesh(decsg(gd, sf, ns), 'Hmax', hmax);

g = Grid(p, e, t);

end
