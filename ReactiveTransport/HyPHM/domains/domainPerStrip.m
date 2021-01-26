%> @file domainPerStrip.m Generation of three domains.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @code [OMe, OM, Y] = domainPerStrip(a, b, phi, e, L, hmax) @endcode
%>
%> This script returns three domains: @c OMe, @c OM and @c Y.  The domain @c OMe
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
%>  @c OM , ie, @f$\Omega@f$ has the same geometry without the inner heterogeneities.  The unit
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
%> @image html OMe.jpg
%> @f$\Omega_\varepsilon@f$
%> @image html OM.jpg
%> @f$\Omega@f$
%> @image html Y.jpg
%> @f$Y@f$
%>
%> @sa FoldedGrid, Grid


%> @param  a          length of axis a (in unit cell) [scalar]
%> @param  b          length of axis b (in unit cell) [scalar]
%> @param  phi        time step length [scalar]
%> @param  e          epsilon [scalar]
%> @param  L          length of domain, has to be a multiple of e [scalar]
%> @param  hmax       mesh fineness of the Grid Y, others will be scaled [scalar]
%>
%> @retval  OMe     semi-periodic grid in HyPHM format with heterogeneities [FoldedGrid]
%> @retval  OM       semi-periodic grid in HyPHM format [FoldedGrid]
%> @retval  Y       periodic grid in HyPHM format [FoldedGrid]
function [OMe, OM, Y] = domainPerStrip(a, b, phi, e, L, hmax)

assert(isscalar(a))
assert(isscalar(b))
assert(isscalar(phi))
assert(isscalar(e))
assert(isscalar(L))
assert(a < 1)
assert(b < 1)
assert(rem(L, e) == 0, 'HyPHM: L has to be a multiple of e.')

numhet = L / e; % number of horizontal heterogeneities

%%%%%%%%%%%%%%%

%% Unit Cell %%
%%%%%%%%%%%%%%%

gd = [3.0000, 4.0000; ... % geometry description
    4.0000, 0.5000; ...
    0, 0.5000; ...
    1.0000, a / 2; ...
    1.0000, b / 2; ...
    0, phi; ...
    0, 0; ...
    0, 0; ...
    1.0000, 0; ...
    1.0000, 0];

sf = 'square-ellipse'; % set formula
ns = names2ns('square', 'ellipse'); % name space
[p, ed, t] = initmesh(decsg(gd, sf, ns), 'Hmax', hmax);
g = Grid(p, ed, t);
Y = FoldedGrid(g);

%%%%%%%%%%%%%%%%%

%% Homo Domain %%
%%%%%%%%%%%%%%%%%

gd = [3; 4; 0; L; L; 0; 0; 0; e; e]; % geometry description
sf = 'rect'; % set formula
ns = names2ns('rect'); % name space
[p, ed, t] = initmesh(decsg(gd, sf, ns), 'Hmax', e*hmax);
g = Grid(p, ed, t);
OM = FoldedGrid(g, 1); % only vertically folded

%%%%%%%%%%%%%%%%%%%

%% Hetero Domain %%
%%%%%%%%%%%%%%%%%%%

gd = [3; 4; 0; L; L; 0; 0; 0; e; e]; % geometry description
sf = 'rect'; % set formula
nslist = {'rect'};

for kL = 1:numhet
    gd = horzcat(gd, [4; (kL - .5) * e; .5 * e; e * a / 2; e * b / 2; phi; 0; 0; 0; 0]); %#ok<AGROW>
    sf = horzcat(sf, ['-', 'cell', int2str(kL)]); %#ok<AGROW>
    nslist{kL+1} = ['cell', int2str(kL)];
end
ns = names2ns(nslist);
[p, ed, t] = initmesh(decsg(gd, sf, ns), 'Hmax', e*hmax);
g = Grid(p, ed, t);
OMe = FoldedGrid(g, 1);

end
