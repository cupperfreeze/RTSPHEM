%>  @file +P1BubbleP1Bubble/visslice.m Visualizes a slice for the type @c P1Bubble.

function visslice(grid, name, filename, data)

% if grid is folded, re-expand to rectangular grid
if isa(grid, 'FoldedGrid')
    V0T = grid.rectGrid.V0T;
    coordV = grid.rectGrid.coordV;
    numV = grid.rectGrid.numV;
    if size(data, 1) == grid.numV + grid.numT
        data = data(grid.mapV, :);
    end % else do nothing, since Tk -> Tk when folded
else
    V0T = grid.V0T;
    coordV = grid.coordV;
    numV = grid.numV;
end
data = full(data(1:numV, :));
vtkquiver(coordV(:, 1), coordV(:, 2), zeros(numV, 1), ...
    data(:, 1), data(:, 2), zeros(numV, 1), name, filename);

end
