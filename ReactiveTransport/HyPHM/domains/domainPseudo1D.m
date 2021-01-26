%> @file domainPseudo1D.m Mesh for onedimensional simulation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> [Grid] = domainPseudo1D(N)
%> Rectangular bounded domain with lower left corner (0, 0), upper
%> right corner (1, 1/N), and mesh fineness 1/N.
%>
%>  @param  N  Mesh fineness in x direction [scalar]
%>  @retval g  The grid in HyPHM format [Grid]
%>
%> See also HyPHM, AbstractGrid, FoldedGrid, Grid

function g = domainPseudo1D(N)

assert(isscalar(N) && N > 0)

coordV = [[(0:1 / N:1)'; (0:1 / N:1)'], [zeros(N + 1, 1); ones(N + 1, 1) / N]];
V0T = [[1:N, 1:N]', [2:N + 1, N + 3:2 * N + 2]', [N + 3:2 * N + 2, N + 2:2 * N + 1]'];

g = Grid(coordV, V0T);

g.idE(g.baryE(:, 2) == 0) = 1; % bottom edges
g.idE(g.baryE(:, 2) == 1/N) = 1; % top edges
g.idE(g.baryE(:, 1) == 0) = 3; % left edge
g.idE(g.baryE(:, 1) == 1) = 2; % right edge
printline(1, 'Boundary IDs set!')

end
