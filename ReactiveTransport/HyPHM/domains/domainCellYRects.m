%> @file domainCellY_rects.m Generation of a unit square-bounded domain with inner rectangular obstacles.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @brief Generation of a unit square-bounded domain with inner rectangular obstacles.
%>
%> @code [Grid] = domainCellY_rects(a1, b1, c1, phi1, a2, b2, c2, phi2, ..., hmax) @endcode
%>
%> Unit square-bounded domain @f$[0,1]\times[0,1]@f$ with inner rectangular obstacles.  A single rectangle has
%> sides of lengths a and b and is rotated by phi around the center c.  The grid is not yet
%> folded.
%> @n @n
%> If You want to use periodic boundaries, do
%> @code
%>   fg = FoldedGrid(domainCellY_rect(...))
%> @endcode
%>
%> See also HyPHM, AbstractGrid, FoldedGrid, Grid
%>
%> @param    ak     Length of side a of the kth rectangle [scalar].
%> @param    bk     Length of side b of the kth rectangle [scalar].
%> @param    ck     Center of the kth rectangle [2 x 1].
%> @param    phik   Rotation angle of of the kth rectangle [scalar].
%> @param    hmax   Fineness of the grid, required by mesh-generator [scalar].
%> @retval   g      Instance of the Grid.

function g = domainCellYRects(varargin)

assert(mod(nargin, 4)-1 == 0, 'Wrong arguments.  Please see manual.')
hmax = varargin{end};
numR = (nargin - 1) / 4; % number of interior rectangles (this is an interger due to assert)
namelist = cell(numR+1, 1);
namelist{1} = 'square';

gd = [3.0000; ... % geometry description (exterior boundary)
    4.0000; ...
    0; ...
    1.0000; ...
    1.0000; ...
    0; ...
    0; ...
    0; ...
    1.0000; ...
    1.0000];

sf = 'square'; % set formula


for k = 1:numR

    a = varargin{4*(k - 1)+1};
    b = varargin{4*(k - 1)+2};
    c = varargin{4*(k - 1)+3};
    phi = varargin{4*(k - 1)+4};

    assert(isscalar(a))
    assert(isscalar(b))
    assert(isvector(c) && length(c) == 2)
    assert(isscalar(phi))
    assert(a < 1 && a > 0)
    assert(b < 1 && b > 0)


    Q = [cos(phi), -sin(phi); sin(phi), cos(phi)]; % rotation matrix
    vec1 = [a / 2; 0]; % vector from center to barycenter of right side
    vec2 = [0; b / 2]; % vector from center to barycenter of top side
    vec1 = Q * vec1; % rotation by phi
    vec2 = Q * vec2; % rotation by phi
    cornBL = c - vec1 - vec2;
    cornTL = c - vec1 + vec2;
    cornTR = c + vec1 + vec2;
    cornBR = c + vec1 - vec2;

    gd_ext = [3; ... % geometry description
        4; ...
        cornBL(1); ...
        cornTL(1); ...
        cornTR(1); ...
        cornBR(1); ...
        cornBL(2); ...
        cornTL(2); ...
        cornTR(2); ...
        cornBR(2)];
    gd = [gd, gd_ext];

    sf = [sf, '-rect', int2str(k)]; % set formula
    namelist{k+1} = ['rect', int2str(k)]; % name space


end

ns = names2ns(namelist); % name space

[p, e, t] = initmesh(decsg(gd, sf, ns), 'Hmax', hmax);

g = Grid(p, e, t);

end
