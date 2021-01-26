%>  @file +P1P1/visslice.m Visualizes a slice for the type @c P1P1.

function visslice(grid, name, filename, data)

vtkquiver(grid.coordV(:, 1), grid.coordV(:, 2), [], ...
    data(:, 1), data(:, 2), [], name, filename);

end
