%>  @file +P1/visslice.m Visualizes a slice for the type @c P1.

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

% vtktrisurf checks if data is on elements or vertices
vtktrisurf(V0T, coordV(:, 1), coordV(:, 2), zeros(numV, 1), data, name, filename);

end
