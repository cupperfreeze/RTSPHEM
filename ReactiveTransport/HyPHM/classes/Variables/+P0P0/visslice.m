%> @file +P0P0/visslice.m Visualizes a slice for the type @c P0P0.

function visslice(grid, name, filename, data)

vtkquiver(grid.baryT(:, 1), grid.baryT(:, 2), [], ...
    data(:, 1), data(:, 2), [], name, filename);

end
