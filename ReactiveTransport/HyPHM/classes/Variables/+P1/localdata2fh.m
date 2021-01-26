%> @file +P1/localdata2fh.m Returns a function_handle @f$f:T\rightarrow \mathbb{R}@f$
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @b Example (Visualization)
%>
%> We want to visualize a discrete function @f$w\in P_1(\mathfrak{T})@f$ (see Variables.type)
%> with local coordinates @f$(w_1,w_2,w_3)^T@f$ on the triangle @f$T_1@f$ and
%> zeros elsewhere.  The returned function_handle can be sampled by meshgrid everywhere,
%> since zero is returned in @f$\Omega\setminus T_1@f$.
%> @code
%> g = Grid([0,0;1,0;1,1;0,1], [1,2,3;1,3,4]);
%> kT = 1;             % triangle no 1
%> coordsT = [3;1;-1]; % coordinates
%> fun = P1.localdata2fh(g, coordsT, 1)
%> [X,Y] = meshgrid(0:0.02:1);
%> for k = 1 : length(X)
%>   for ell = 1 : length(X)
%>     Z(k, ell) = fun([X(k, ell); Y(k, ell)]);
%>   end
%> end
%> surf(X, Y, Z)
%> @endcode
%>
%> @image html localP1function.jpg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @param grid Instance of Grid.
%> @param coordsT Coordinates of @f$w_h\in P_1@f$ [3 x 1].
%> @param kT Number of considered triangle @f$T@f$.

function fun = localdata2fh(grid, coordsT, kT)

assert(length(coordsT) == 3, 'HyPHM: Three coordinates (for each vertex) required [3x1].')

fun = @(X) dot(coordsT, P1.localbasis(grid, kT, X));

end
