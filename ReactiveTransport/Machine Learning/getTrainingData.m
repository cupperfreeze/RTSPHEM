% obtain numberOut x  innerIterations image-label pairs
% images are of size [AimSize x AimSize x 1] and extracted from image

function [outputIm, outputData, physics] = getTrainingData(numberOut, AimSize, innerIterations, image)

outputIm = zeros(AimSize, AimSize, 1, numberOut*innerIterations, 'int8'); % Image storage
outputData = zeros(numberOut*innerIterations, 4); % Label storage
physics = zeros(numberOut*innerIterations, 2); % physical parameter storage
TotalSize = size(image, 1);

cellGrid = FoldedCartesianGrid(2, [-0.5, 0.5, -0.5, 0.5], (AimSize - 1)*[1, 1]);
GridHyPHM = Grid(cellGrid.coordinates, cellGrid.triangles);

for i = 1:numberOut
    % randomly extract sample
    coord = AimSize / 2 + ceil(rand(2, 1)*(TotalSize - AimSize));
    layer = ceil(rand(1, 1)*TotalSize);
    currentImage = double(image((coord(1)-AimSize/2+1):(coord(1) + AimSize / 2), (coord(2) - AimSize / 2 + 1):(coord(2) + AimSize / 2), layer));
    %imagesc(currentImage)
    % convert to levelSet representation
    currentImage(currentImage <= 128) = -1;
    currentImage(currentImage > 128) = 1;
    levelSet = reinitializeLevelSet(cellGrid, -currentImage(:), true, inf);
    if mod(i, 50) == 0
        i
    end
    for j = 1:innerIterations
        numSteps = floor(1+6*rand(1));
        % Perform level-Set evolution
        for k = 1:numSteps
            levelSet = levelSetEquationTimeStep(1/AimSize, 0, levelSet, ...
                cellGrid, -0.5*(j > 1));
        end
        %figure
        %contour(reshape(levelSet,AimSize,AimSize),[0,1])
        % calculate permeability labels
        %

        permeability = computePermeabilityTensor(GridHyPHM, ...
            levelSet, 'Bubble');
        %          [ ~, ~, ~, triangleVolumes, triangleSurfaces ] ...
        %         = assembleCellProblem( cellGrid, levelSet );

        % physics((i-1)*innerIterations+j,:) = [sum(triangleVolumes), sum(triangleSurfaces)];
        outputIm(:, :, 1, (i - 1)*innerIterations+j) = reshape(levelSet > 0, AimSize, AimSize);
        outputData((i-1)*innerIterations+j, :) = permeability(:);
    end


end
