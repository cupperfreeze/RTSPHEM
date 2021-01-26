%> @file +RT0/localdata2fh.m Returns a function_handle @f$\vec{f}:T\rightarrow \mathbb{R}^2@f$
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @b Example (Visualization)
%>
%> @code
%> g = Grid([0,0;1,0;1,1;0,1], [1,2,3;1,3,4]); % unit square with two triangles
%> dx = 25E-3; % sample fineness
%> [X, Y] = meshgrid(0:dx:1);
%>
%> dim = length(X);
%>
%> % global solution: -2, -2.5, -3, -1, 1, cf g.visualize('numE')
%>
%> hold on
%>
%> % triangle T1
%> U = zeros(dim);  V = zeros(dim);
%> f = RT0.localdata2fh(g, [-3;-2.5;-2], 1); % local solutions as function_handle
%> for k = 1 : 1/dx+1
%>   for ell = 1 : 1/dx+1
%>     tmp = f([X(k, ell); Y(k, ell)]);
%>     U(k, ell) = tmp(1);
%>     V(k, ell) = tmp(2);
%>   end
%> end
%> quiver(X,Y,U,V)
%>
%> % triangle T2
%> f = RT0.localdata2fh(g, [1;-1;-2.5], 2); % local solutions as function_handle
%> for k = 1 : 1/dx+1
%>   for ell = 1 : 1/dx+1
%>     tmp = f([X(k, ell); Y(k, ell)]);
%>     U(k, ell) = tmp(1);
%>     V(k, ell) = tmp(2);
%>   end
%> end
%> quiver(X,Y,U,V)
%> @endcode
%>
%  @image html localP2function.jpg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @param grid Instance of Grid.
%> @param coordsT Coordinates of @f$\vec{w}_h\in \vec{RT}_0@f$ [3 x 1].
%> @param kT Number of considered triangle @f$T@f$.

function fun = localdata2fh(grid, coordsT, kT)

assert(isequal(size(coordsT), [3, 1]), 'HyPHM: Three coordinates (for each edge) required [3x1].')

fun = @(X) (coordsT' * RT0.localbasis(grid, kT, X))';

end
