%> @file +P2/localdata2fh.m Returns a function_handle @f$f:T\rightarrow \mathbb{R}@f$
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @b Example (Visualization)
%>
%> See P1/localdata2fh.m.
%>
%> @image html localP2function.jpg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @param grid Instance of Grid.
%> @param coordsT Coordinates of @f$w_h\in P_2@f$ [6 x 1].
%> @param kT Number of considered triangle @f$T@f$.

function fun = localdata2fh(grid, coordsT, kT)

assert(length(coordsT) == 6, 'HyPHM: Six coordinates (for each vertex and edge) required [6x1].')

fun = @(X) dot(coordsT, P2.localbasis(grid, kT, X));

end
