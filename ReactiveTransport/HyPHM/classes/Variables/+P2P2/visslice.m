%> @file +P2P2/visslice.m Visualizes a slice for the type @c P2P2.

function visslice(grid, name, filename, data)

data = full(data);
vtkquiver([grid.coordV(:, 1); grid.baryE(:, 1)], [grid.coordV(:, 2); grid.baryE(:, 2)], [], ...
    data(:, 1), data(:, 2), [], name, filename);

end
