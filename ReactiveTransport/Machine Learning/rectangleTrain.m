% obtain numberOut many image-label pairs of artificial geometries
% consisting of rectangles and ellipses

function [outputIm, outputData, physics] = rectangleTrain(numberOut, AimSize)

outputIm = zeros(AimSize, AimSize, 1, numberOut, 'int8');
outputData = zeros(numberOut, 4);
physics = zeros(numberOut, 2);

cellGrid = FoldedCartesianGrid(2, [-0.5, 0.5, -0.5, 0.5], (AimSize - 1)*[1, 1]);
coord = cellGrid.coordinates;
GridHyPHM = Grid(cellGrid.coordinates, cellGrid.triangles);

for i = 1:numberOut
    image = -ones(AimSize, AimSize);

    % add random rectangles and ellipses to empty image
    for j = 1:ceil(3*rand(1))
        image = rectangle(image, coord, -0.4+0.8*rand(1, 2), 0.025+0.3*rand(1), 0.025+0.3*rand(1), AimSize);
    end
    for j = 1:ceil(2*rand(1))
        image = ellipse(image, coord, -0.4+0.8*rand(1, 2), 0.05+0.2*rand(1), 0.05+0.2*rand(1), AimSize);
    end
    %imagesc(image)

    % convert to level-set representation and calculate permeability
    % labels
    levelSet = reinitializeLevelSet(cellGrid, image(:), true, inf);
    permeability = computePermeabilityTensor(GridHyPHM, ...
        levelSet, 'Bubble');
    [~, ~, ~, triangleVolumes, triangleSurfaces] ...
        = assembleCellProblem(cellGrid, levelSet);
    physics(i, :) = [sum(triangleVolumes), sum(triangleSurfaces)];
    outputIm(:, :, 1, i) = reshape(levelSet > 0, AimSize, AimSize);
    outputData(i, :) = permeability(:);

end

end

% add rectangle
function out = rectangle(in, coord, mid, x, y, AimSize)
new = -ones(size(in));
new(abs(coord(:, 1) - mid(1)) <= x & abs(coord(:, 2) - mid(2)) <= y) = 1;
out = reshape(max([in(:), new(:)], [], 2), AimSize, AimSize);
end

% add ellipse
function out = ellipse(in, coord, mid, x, y, AimSize)
new = -ones(size(in));
new(((coord(:, 1) - mid(1)).^2./x.^2+(coord(:, 2) - mid(2)).^2./y.^2) < 1) = 1;
out = reshape(max([in(:), new(:)], [], 2), AimSize, AimSize);
end