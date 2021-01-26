%>  @file +CR/visslice.m Visualizes a slice for the type @c CR [#E x 1].

function visslice(grid, name, filename, data)


% if grid is folded, re-expand to rectangular grid
if isa(grid, 'FoldedGrid')
    V0T = grid.rectGrid.V0T;
    coordV = grid.rectGrid.coordV;
    numV = grid.rectGrid.numV;
    if length(data) == grid.numV
        data = data(grid.mapV);
    end % else do nothing, since Tk -> Tk when folded
else
    V0T = grid.V0T;
    coordV = grid.coordV;
    numV = grid.numV;
end

% here we must extrapolate the values on edges to the vertices.  Let v_k
% denote the vertices of a fixed triangle T and let e_k denote the
% baricenters of the edges opposite to v_k.  Let m_k denote the midpoint
% between v_k and e_k.  The point is too a midpoint of e_l and e_m,
% l~=m~=k.  Not we can use the extrapolation rule
% f(v_k) = 2*f(m_k) - f(e_k) = f(e_m) + f(e_l) - f(e_k).

% @todo This needs to be vectorized!
visdata = zeros(grid.numT, 3);
for kT = 1:grid.numT
    for k = 1:3
        l = mod(k, 3) + 1;
        m = mod(k+1, 3) + 1;

        visdata(kT, k) = data(grid.E0T(kT, l)) + data(grid.E0T(kT, m)) - data(grid.E0T(kT, k));
    end
end

% vtktrisurf checks if data is on elements or vertices
vtktrisurf(V0T, coordV(:, 1), coordV(:, 2), zeros(numV, 1), visdata, name, filename);

end
