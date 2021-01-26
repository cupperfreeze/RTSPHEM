% [Grid] = domainCellYCraggy(diam, cragfac, hmax)
%
% Unit square-bounded domain with inner cragged hole.  The hole has
% an average diameter of diam.  The grid is not yet folded.  If You want
% to use periodic boundaries, do
% fg = FoldedGrid(domainCellYCraggy(diam, cragfac, hmax)).
%
% Input
%   diam       [ scalar ]  diameter of inner hole
%   cragfac    [ scalar ]  amount of craggyfication, zero is none
%   hmax       [ scalar ]  fineness of grid
%
% Output
%   g          [ Grid ]    the grid in HyPHM format
%
% See also HyPHM, AbstractGrid, FoldedGrid, Grid, domainCellY
%
% Copyright 2009, 2010 F. Frank, Chair of Applied Mathematics I,
%           Department of Mathematics, University of Erlangen-Nuremberg,
%           G E R M A N Y

function g = domainCellYCraggy(diam, cragfac, hmax)

assert(isscalar(diam))
assert(isscalar(cragfac))
assert(diam > 0 && diam < 1)
assert(cragfac >= 0)

theta = 0:.1:2 * pi; % rotation angle
numS = length(theta); % number of samples
radius = diam / 2;
craggyradius = (1 + cragfac * (rand(numS, 1) - .5)) * radius;

gd_rect = [3; 4; 0; 1; 1; 0; 0; 0; 1; 1];
gd_grain = [2; numS; ...
    .5 + craggyradius .* sin(theta)'; ...
    .5 + craggyradius .* cos(theta)'];

gd = [[gd_rect; zeros(2 * numS + 2 - 10, 1)], gd_grain]; % geometry description

sf = 'square-grain'; % set formula
ns = names2ns('square', 'grain'); % name space

[p, e, t] = initmesh(decsg(gd, sf, ns), 'Hmax', hmax);

g = Grid(p, e, t);

printline(-1, 'Note that the boundary is generated RANDOMLY each time the grid is generated!')

end
