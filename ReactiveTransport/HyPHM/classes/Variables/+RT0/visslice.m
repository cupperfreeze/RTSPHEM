%> @file +RT0/visslice.m Visualizes a slice for the type @c RT0.

function visslice(grid, name, filename, data)

flux = RT0.RT0toP1P1slice(grid, data);

switch isa(grid, 'FoldedGrid')

    case true % loop over each single triangle to ged periodic (unfolded) data

        udata = zeros(grid.numT, 3); % data on each vertex for each triangle
        vdata = zeros(grid.numT, 3);

        for kT = 1:grid.numT
            udata(kT, 1:3) = flux(grid.V0T(kT, :), 1);
            vdata(kT, 1:3) = flux(grid.V0T(kT, :), 2);
        end

        verts1 = squeeze(grid.coordV0T(:, :, 1));
        verts2 = squeeze(grid.coordV0T(:, :, 2));

        vtkquiver(verts1(:), verts2(:), [], ...
            udata(:), vdata(:), [], name, filename);


    case false % use global coordinates

        vtkquiver(grid.coordV(:, 1), grid.coordV(:, 2), [], ...
            flux(:, 1), flux(:, 2), [], name, filename);

end

end
