%> @file +P2P2/localdata2fh.m Returns a function_handle @f$f:T\rightarrow \mathbb{R}^2@f$
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> Depends on P2/localbasis.m.
%>
%> @b Example (Visualization)
%>
%> See P1/localdata2fh.m.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @param grid Instance of Grid.
%> @param coordsT Coordinates of @f$w_h\in P_1@f$ @f$[6 \times 2]@f$.
%> @param kT Number of considered triangle @f$T@f$.

function fun = localdata2fh(grid, coordsT, kT)

assert(isequal(size(coordsT), [6, 2]), 'HyPHM: Six coordinates (for each vertex and edge) per component required [6x2].')

fun = @(X) [dot(coordsT(:, 1), P2.localbasis(grid, kT, X)); ...
    dot(coordsT(:, 2), P2.localbasis(grid, kT, X))];

end
