%> @file domainPerStrip3.m Generation of three domains.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @brief Generation of three domains.
%>
%> @code [ghetero, ghomo, gcell] = domainPerStrip(a, b, phi, e, L, hmax) @endcode
%>
%> This script returns three domains: @c ghetero, @c ghomo and @c gcell.  The domain @c ghetero
%> representing @f$\Omega_\varepsilon@f$ has the geometry
%>
%> @code
%>   perperperperperperperperperperperperperperperper   e
%>   +                                              +   |
%>   +   ++      ++                           ++    +   |
%>   +  +  +    +  +                         +  +   +   |
%>   +  +   +   +   +                        +   +  +   |
%>   +   +  +    +  +                         +  +  +   |
%>   +    ++      ++                           ++   +   |
%>   +                                              +   |
%>   perperperperperperperperperperperperperperperper   0
%>
%>   0----------------------------------------------L
%>   0-------e
%>  @endcode
%>
%>  @c ghomo , ie, @f$\Omega@f$ has the same geometry without the inner heterogeneities.  The unit
%>  cell @f$Y@f$  looks like
%>
%>  @code
%>   perperperperperperper   1
%>   per               per   |
%>   per      ++       per   |
%>   per     +  +      per   |
%>   per     +   +     per   |
%>   per      +  +     per   |
%>   per       ++      per   |
%>   per               per   |
%>   perperperperperperper   0
%>
%>   0-------------------1
%> @endcode
%>
%> The heterogeneity itself is an ellipse with axes length @c a and @c b, rotated
%>  by @c phi.
%>
%> To check the IDs, use <tt>FoldedGrid.visualize('idE')</tt>.
%>
%> @image html ghetero.jpg
%> @f$\Omega_\varepsilon@f$
%> @image html ghomo.jpg
%> @f$\Omega@f$
%> @image html gcell.jpg
%> @f$Y@f$
%>
%> @sa FoldedGrid, Grid


%> @param  e          epsilon [scalar]: scaling factor
%> @param  L          length of domain, has to be a multiple of e [scalar]
%> @param  hmax       mesh fineness of the Grid gcell, others will be scaled [scalar]
%>
%> @retval  ghetero     semi-periodic grid in HyPHM format with heterogeneities [FoldedGrid]
%> @retval  ghomo       semi-periodic grid in HyPHM format [FoldedGrid]
%> @retval  gcell       periodic grid in HyPHM format [FoldedGrid]
function [ghetero, ghomo, gcell] = domainPerStrip3(e, L, hmax)

% mxk         x coordinate of the midpoint of ellipse k (in unit cell) [scalar]
% myk         y coordinate of the midpoint of ellipse k (in unit cell) [scalar]
% ak          first radius of ellipse k (in unit cell) [scalar]
% bk          second radius of ellipse k (in unit cell) [scalar]
% phik        rotation angle of ellipse k (in unit cell) [scalar]

% Data for Ellipse 1
mx1 = 0.25;
my1 = 2 / 3;
a1 = 0.5;
b1 = 0.25;
phi1 = pi / 4;

% Data for Ellipse 2
mx2 = 0.5;
my2 = 0.25;
a2 = 0.5;
b2 = 0.25;
phi2 = -pi / 4;

% Data for Ellipse 3
mx3 = 0.75;
my3 = 2 / 3;
a3 = 0.5;
b3 = 0.25;
phi3 = pi / 2;

assert(isscalar(a1))
assert(isscalar(b1))
assert(isscalar(a2))
assert(isscalar(b2))
assert(isscalar(a3))
assert(isscalar(b3))
assert(isscalar(phi1))
assert(isscalar(phi2))
assert(isscalar(phi3))
assert(isscalar(e))
assert(isscalar(L))
assert(a1 < 1)
assert(b1 < 1)
assert(a2 < 1)
assert(b2 < 1)
assert(a3 < 1)
assert(b3 < 1)
assert(rem(L, e) == 0, 'HyPHM: L has to be a multiple of e.')

numhet = L / e; % number of horizontal heterogeneities

%%%%%%%%%%%%%%%

%% Unit Cell %%
%%%%%%%%%%%%%%%

gd = [3.0000, 4.0000, 4.0, 4.0; ... % geometry description
    4.0000, mx1, mx2, mx3; ...
    0, my1, my2, my3; ...
    1.0000, a1 / 2, a2 / 2, a3 / 2; ...
    1.0000, b1 / 2, b2 / 2, b3 / 2; ...
    0, phi1, phi2, phi3; ...
    0, 0, 0, 0; ...
    0, 0, 0, 0; ...
    1.0000, 0, 0, 0; ...
    1.0000, 0, 0, 0];

ns = names2ns('square', 'ellipse1', 'ellipse2', 'ellipse3'); % name space
sf = 'square-ellipse1-ellipse2-ellipse3'; % set formula
[p, ed, t] = initmesh(decsg(gd, sf, ns), 'Hmax', hmax);
g = Grid(p, ed, t);
gcell = FoldedGrid(g);

%%%%%%%%%%%%%%%%%

%% Homo Domain %%
%%%%%%%%%%%%%%%%%

gd = [3; 4; 0; L; L; 0; 0; 0; e; e]; % geometry description
sf = 'rect'; % set formula
ns = names2ns('rect'); % name space
[p, ed, t] = initmesh(decsg(gd, sf, ns), 'Hmax', e*hmax);
g = Grid(p, ed, t);
ghomo = FoldedGrid(g, 1); % only vertically folded

% %%%%%%%%%%%%%%%%%%%

%% Hetero Domain %%
%%%%%%%%%%%%%%%%%%%

gd = [3; 4; 0; L; L; 0; 0; 0; e; e]; % geometry description
nslist = {'rect'};
sf = 'rect'; % set formula


for k = 1:numhet
    gd = horzcat(gd, [4; (k - (1 - mx1)) * e; my1 * e; e * a1 / 2; e * b1 / 2; phi1; 0; 0; 0; 0], ...
        [4; (k - mx2) * e; my2 * e; e * a2 / 2; e * b2 / 2; phi2; 0; 0; 0; 0], ...
        [4; (k - (1 - mx3)) * e; my3 * e; e * a3 / 2; e * b3 / 2; phi3; 0; 0; 0; 0]); %#ok<AGROW>
    sf = horzcat(sf, ['-', 'cellA', int2str(k)], ['-', 'cellB', int2str(k)], ['-', 'cellC', int2str(k)]); %#ok<AGROW>
    nslist{3*k-1} = ['cellA', int2str(k)];
    nslist{3*k} = ['cellB', int2str(k)];
    nslist{3*k+1} = ['cellC', int2str(k)];
end
ns = names2ns(nslist);

[p, ed, t] = initmesh(decsg(gd, sf, ns), 'Hmax', e*hmax);
g = Grid(p, ed, t);
ghetero = FoldedGrid(g, 1);

end
