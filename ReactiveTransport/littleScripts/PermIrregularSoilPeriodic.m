% Evolve geometry according to two-phase solid with different normal
% interface velocities. Porosity and permeability of the arising unit cell
% geometries are calculated
close all
global EPS
EPS = eps;

dimension = 2;
numPartitions = 64; % Number of partitions in each direction
endTime = 1.2;
timeStepSize = 1 / numPartitions;

cellGrid = FoldedCartesianGrid(dimension, kron(ones(1, dimension), [-0.5, 0.5]), numPartitions*ones(1, dimension));
coord = cellGrid.coordinates;

initialLevelSetFunc = @(x) abs(x(1)) - 0.1;
% lsf = @(x) min( abs( x(1) + 0.25 ) - 0.05, abs( x(1) - 0.25 ) - 0.05 );
coordCell = mat2cell(coord, ones(1, cellGrid.nodes), dimension);
initialData = cellfun(initialLevelSetFunc, coordCell);


timeCell = num2cell(ones(cellGrid.nodes, 1));
structure = cellfun(@speedFunction, timeCell, coordCell);
figure
imshow(reshape(structure, numPartitions + 1, numPartitions + 1)', 'InitialMagnification', 'fit')
clear coordCell;
X = cellGrid.reshape(coord(:, 1));
Y = cellGrid.reshape(coord(:, 2));
clear coord;

plotInitial = cellGrid.reshape(initialData);

% speed = @(t,x) -35 * (x(1)-0.3).^2 * ( x(1) + 0.5 ).^2 * ( x(1) <= 0.3 );
% variableSpeed = @(t,x) speed(t,x) * ( ( t <= 0.5*endTime ) - ( t > 0.5*endTime ) );

fig = figure;
[~, defHandle] = contour(X, Y, plotInitial, zeros(1, 2), 'LineColor', 'r', ...
    'LineWidth', 3);
set(gca, 'DataAspectRatio', [1, 1, 1]);
% xlabel('x');
% ylabel('y');
handle = copyobj(defHandle, gca);
set(handle, 'Zdata', plotInitial, 'LineColor', 'b', 'LineWidth', 3);
set(gca, 'xtick', []);
set(gca, 'ytick', []);


toPrint = false;
printFileName = ['PermIrregularSoil', 'TimeStep_'];
if (toPrint)
    print(fig, '-depsc', [printFileName, num2str(0), '.eps']);
end
toSave = false; % Flag: save variables to disk
speedFun = @speedFunction;


levelSet = solveLevelSetEquationOld(cellGrid, initialData, @speedFunction, ...
    [], endTime, timeStepSize, 'ContourHandle', handle, 'SaveTime', 'all', ...
    'ContourFigure', fig, 'Print', toPrint, 'PrintFileName', printFileName, ...
    'Reinitialize', false);

clear velocity;
clear fig handle defHandle;

numTimeSlices = size(levelSet, 2);
permeabilityTensors = NaN(4, numTimeSlices);
porosity = NaN(numTimeSlices, 1);
coord = cellGrid.coordinates;
helpGridHyPHM = Grid(coord, cellGrid.triangles);


for ts = 1:numTimeSlices
    disp(['Time step ', num2str(ts)]);
    [cellProblemSystemMatrix, rhs, isDoF, triangleVolumes, triangleSurfaces] ...
        = assembleCellProblem(cellGrid, levelSet(:, ts));
    porosityValues(ts) = sum(triangleVolumes);
    surfaceValues(ts) = sum(triangleSurfaces);
    put = computePermeabilityTensor(helpGridHyPHM, levelSet(:, ts));
    permeabilityTensors(:, ts) = put(:);
end
clear ts diffusion porosity;

% copyfile( filePath, [ folderName, 'Data' ] );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function normalVelocity = speedFunctionTemplate(t, x)

calciteSpeed = 0.5;
dolomiteSpeed = 0.05;

if (any(abs(x) > 0.5))
    normalVelocity = 0;
    return;
end

if (x(1) < -0.15 && x(2) > -0.5 * x(1) + 0.475)
    normalVelocity = dolomiteSpeed;
    return;
end

if abs(x(1)) > 0.3 & x(2) < 0
    normalVelocity = dolomiteSpeed;
    return;
end

if (x(1) < 0 && ((x(2) < 0.3 * x(1) + 0.4) && (x(2) > 0.2 * x(1) + 0.3)))
    normalVelocity = dolomiteSpeed;
    return;
end

if ((x(2) > 0.3 * x(1) + 0.0) && (x(2) < 0.3 * x(1) + 0.25))
    normalVelocity = dolomiteSpeed;
    return;
end

if (x(2) < (0.3 * (x(1) < 0) + 0.05 * (x(1) >= 0)) * x(1) - 0.1)
    if (x(2) > (-0.3 * (x(1) < 0) + -0.25 * (x(1) >= 0)))
        normalVelocity = dolomiteSpeed;
        return;
    end
end

if (x(2) < (0.2 * (x(1) < 0)) * x(1) - 0.35)
    normalVelocity = dolomiteSpeed;
    return;
end

normalVelocity = calciteSpeed;

end

function normalVelocity = speedFunction(t, x)

y = x;
y(x >= 0) = 2 * y(x >= 0) - 0.5;
y(x < 0) = -2 * y(x < 0) - 0.5;

normalVelocity = speedFunctionTemplate(t, y);

end
