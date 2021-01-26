%>  @file +P1disc/visslice.m Visualizes a slice for the type @c P1disc.

function visslice(grid, name, filename, data)

% data [3*#T, 1]  -->  [#T, 3]
data = reshape(data, 3, grid.numT)';
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

% vtktrisurf checks if data is on elements or vertices
vtktrisurf(V0T, coordV(:, 1), coordV(:, 2), zeros(numV, 1), data, name, filename);

end
