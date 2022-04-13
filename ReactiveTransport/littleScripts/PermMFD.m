% calculate permeability of a microfluidic device
% cf. [4] Section 4.2

global EPS

scale = 30; % number of solid inclusions per dimension
dimension = 2;
numPartitions = 500; % grid resolution
cellGrid = FoldedCartesianGrid(dimension, ...
    kron(ones(1, dimension), [-0.5, 0.5]), ...
    numPartitions*ones(1, dimension));
helpGridHyPHM = Grid(cellGrid.coordinates, cellGrid.triangles);
coord = cellGrid.coordinates;

%% generate level-set function
g = @(x) -1;
initialLevelSetFunc = @(x) g(x);
coordCell = mat2cell(coord, ones(1, cellGrid.nodes), dimension);
levelSet = cellfun(initialLevelSetFunc, coordCell);


% stencil for single inclusion
f = @(x, a) 400 / 1100 / scale - norm(x-a, inf);
f2 = @(x, a) 400 / 1000 / scale - norm(x-a, inf);

for i = 1:scale - 1
    i
    for j = 1:scale - 1
        m = -0.5 + 1 / scale * i; % coordinates of current inclusion
        n = -0.5 + 1 / scale * j;
        if n > 0
            n = n + ((2500 - 300) / 1100) / scale;
            indizes = logical(indizesF([m, n], coord, scale)); % indizes of pixels neighboring this inclusion
            sub = mat2cell(coord(indizes, :), ones(1, sum(indizes)), dimension);
            levelSet(indizes) = max(cellfun(@(x) f2(x, [m, n]), sub), levelSet(indizes)); % add inclusion to level-set function
            continue
        end

        if abs(-0.3*m-0.2-n) > (2500 - 300) / 1100 / scale / 2
            indizes = logical(indizesF([m, n], coord, scale));
            sub = mat2cell(coord(indizes, :), ones(1, sum(indizes)), dimension);
            levelSet(indizes) = max(cellfun(@(x) f(x, [m, n]), sub), levelSet(indizes));
        end

    end
end

% cut-off top/bottom
levelSet = max(cellfun(@(x) 1 - norm([x(1), 1000 / 212 * (x(2) - 0.5)], inf), coordCell), levelSet);
levelSet = max(cellfun(@(x) 1 - norm([x(1), 1000 / 212 * (x(2) + 0.5)], inf), coordCell), levelSet);


permeability = computePermeabilityTensor(helpGridHyPHM, levelSet)

% illustrate final geometry
[a, b] = meshgrid(-0.5:1/numPartitions:0.5, -0.5:1/numPartitions:0.5);
scatter(b(levelSet > 0), a(levelSet > 0), 'filled', 'MarkerFaceColor', 'b')
axis([-0.5, 0.5, -0.5, 0.5])
axis equal
hold on
contour(a, b, reshape(levelSet, numPartitions + 1, numPartitions + 1)', [0, 0], 'linewidth', 1, 'color', 'b')
set(gca, 'xtick', []);
set(gca, 'ytick', []);


function out = indizesF(in, coords, scale)
out = logical(abs(coords(:, 1) - in(1)) < 3/scale & abs(coords(:, 2) - in(2)) < 3/scale);
end

% (80,400)
%   1.0e-07 *
%
%     0.1547   -0.0005
%    -0.0003    0.0004
%
% (80,800)
%
%    1.0e-07 *
%
%     0.1496   -0.0004
%    -0.0004    0.0002

% permeability =
%
%    1.0e-04 *
%
%     0.3870    0.0000
%    -0.0000    0.0000
