%> @file +RT0/localbasis.m Returns the vector of (values of) local basis functions of @f$\mathbb{RT}_0(T)@f$ evaluated in @f$\vec{x}\in T@f$.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> If @f$\vec{x}\notin T@f$ zero ist returned and a warning thrown. Depends on P1/localbasis.m (test whether @f$\vec{x}\in T@f$).
%>
%> The degrees of freedom lie on edges with following order:
%> @code
%>          V3
%>               E1
%>           E2        V2
%>                 E3
%>             V1
%> @endcode
%>
%> @b Example (Usage):
%> See example in P1/localbasis.m.
%>
%>
%> @b Example (Visualization):
%> For this type of finite elements,
%> @f[\vec{\phi}_E\cdot\vec{\nu}_{E'} = \delta_{EE'}~\text{on}~T@f]
%> has to hold, @f$\delta@f$ denoting the Kronecker-Delta (cf @ref notation) and @f$\vec{\nu}_{E'}@f$ the <i>global</i> unit edge normals (for the <i>local</i> edge normals, @f$\vec{\nu}_{ET'}=\sigma_{ET'}\vec{\nu}_{E'}@f$ holds).  This can be done
%> for some sample points, cf eg P2/localbasis.m
%>
%> The visualization can be performed by
%> @code
%>  g = Grid([0,0;1,0;1,1;0,1], [1,2,3;1,3,4]); % unit square with two triangles
%>  dx = 25E-3; % sample fineness
%>  [X, Y] = meshgrid(0:dx:1);
%>  Z = cell(3, 2); % one component for every (local) basis function
%>  for j = 1 : 3
%>    for m = 1 : 2
%>      Z{j,m} = zeros(1/dx+1);
%>    end
%>  end
%>  for k = 1 : 1/dx+1
%>    for ell = 1 : 1/dx+1
%>      tmp = RT0.localbasis(g, 1, [X(k, ell); Y(k, ell)]);
%>      for j = 1 : 3
%>        for m = 1 : 2
%>          Z{j,m}(k, ell) = tmp(j,m);
%>        end
%>      end
%>    end
%>  end
%>
%>  hold on
%>  for j = 1 : 3
%>  subplot(1,3,j)
%>  quiver(X,Y,Z{j,1},Z{j,2})
%>  end
%> @endcode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @param grid  Instance of Grid.
%> @param kT   Global triangle number of @f$T@f$ [scalar].
%> @param x    Evaluation point [2x1].
%> @retval ret Three @f$\mathbb{RT}_0@f$-basis functions  evaluated at @f$\vec{x}\in T@f$ [3x2].


function ret = localbasis(grid, kT, x)
assert(nargin == 3)
assert(isa(grid, 'AbstractGrid'))
assert(isequal(size(x), [2, 1]), 'HyPHM: Third argument has to be [2x1].')
g = grid;

ret = zeros(3, 2);

lam = +P1.localbasis(grid, kT, x); % [3 x 1] barycentric coordinates of x in T or zero if x not in T
if isequal(lam, [0; 0; 0])
    return % return zeros instead of (negative) coordinates when x is not in T
end

coordV = squeeze(g.coordV0T(kT, :, :));

for ell = 1:3
    ret(ell, :) = g.sigE0T(kT, ell) * g.areaE(g.E0T(kT, ell)) / g.areaT(kT) / 2 ...
        * (x - coordV(ell, :)');
end


end
