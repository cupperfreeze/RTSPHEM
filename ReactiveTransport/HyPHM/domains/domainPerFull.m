%> @file domainPerFull.m Generation of three domains @ref Frank2013.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @code [OMe, OM, Y] = domainPerFull(a, b, phi, e, L, H, hmax) @endcode
%>
%> This script returns three domains: @c OMe, @c OM and @c Y.  The domain @c OMe
%> representing @f$\Omega_\varepsilon@f$ has the geometry
%>
%> @code
%>   ++++++++++++++++++++++++++++++++++++++++++++++++   H
%>   +                                              +   |
%>   +   ++      ++                           ++    +   |
%>   +  +  +    +  +                         +  +   +   |
%>   +  +   +   +   +                        +   +  +   |
%>   +   +  +    +  +                         +  +  +   |
%>   +    ++      ++                           ++   +   |
%>   +                                              +   |
%>   +                                              +   |  e
%>   +   ++      ++                           ++    +   |  |
%>   +  +  +    +  +                         +  +   +   |  |
%>   +  +   +   +   +                        +   +  +   |  |
%>   +   +  +    +  +                         +  +  +   |  |
%>   +    ++      ++                           ++   +   |  |
%>   +                                              +   |  |
%>   ++++++++++++++++++++++++++++++++++++++++++++++++   0  0
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
%> @image html domainPerFull-OMe.png
%> @f$\Omega_\varepsilon@f$
%> @image html domainPerFull-OM.png
%> @f$\Omega@f$
%> @image html domainPerFull-Y.png
%> @f$Y@f$
%>
%> @sa FoldedGrid, Grid


%> @param  a          length of axis a (in unit cell) [scalar]
%> @param  b          length of axis b (in unit cell) [scalar]
%> @param  phi        time step length [scalar]
%> @param  e          epsilon [scalar]
%> @param  L          length of domain, has to be a multiple of e [scalar]
%> @param  H          height of domain, has to be a multiple of e [scalar]
%> @param  hmax       mesh fineness of the Grid Y, others will be scaled [scalar]
%>
%> @retval  OMe       grid in HyPHM format with heterogeneities [Grid]
%> @retval  OM        grid in HyPHM format [Grid]
%> @retval  Y       periodic grid in HyPHM format [FoldedGrid]
function [OMe, OM, Y] = domainPerFull(a, b, phi, e, L, H, hmax)

assert(isscalar(a))
assert(isscalar(b))
assert(isscalar(phi))
assert(isscalar(e))
assert(isscalar(L))
assert(isscalar(H))

assert(a < 1)
assert(b < 1)
assert(rem(L, e) == 0, 'HyPHM: L has to be a multiple of e.')
assert(rem(H, e) == 0, 'HyPHM: H has to be a multiple of e.')

numL = L / e; % number of horizontal heterogeneities
numH = H / e; % number of vertical heterogeneities

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

gd = [3; 4; 0; L; L; 0; 0; 0; H; H]; % geometry description
sf = 'rect'; % set formula
ns = names2ns('rect'); % name space
[p, ed, t] = initmesh(decsg(gd, sf, ns), 'Hmax', e*hmax);
OM = Grid(p, ed, t);


%%%%%%%%%%%%%%%%%%%

%% Hetero Domain %%
%%%%%%%%%%%%%%%%%%%

gd = [3; 4; 0; L; L; 0; 0; 0; H; H]; % geometry description
sf = 'rect'; % set formula
nslist = {'rect'};

for kL = 1:numL
    for kH = 1:numH
        gd = horzcat(gd, [4; (kL - .5) * e; (kH - .5) * e; e * a / 2; e * b / 2; phi; 0; 0; 0; 0]); %#ok<AGROW>
        sf = horzcat(sf, ['-', 'cell', int2str(kL), int2str(kH)]); %#ok<AGROW>
        nslist{end+1} = ['cell', int2str(kL), int2str(kH)];
    end
end
ns = names2ns(nslist);
[p, ed, t] = initmesh(decsg(gd, sf, ns), 'Hmax', e*hmax);
OMe = Grid(p, ed, t);


end
