%> @file +C/visslice.m Visualizes a slice for the type @c C, ie a constant on @f$\Omega@f$.

function visslice(grid, name, filename, data)

% vtktrisurf checks if data is on elements or vertices
vtktrisurf(grid.V0T, grid.coordV(:, 1), grid.coordV(:, 2), ...
    zeros(grid.numV, 1), data*ones(grid.numV, 1), name, filename);

end
