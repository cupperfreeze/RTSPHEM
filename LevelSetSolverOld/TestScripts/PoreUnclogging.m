%Level-set simulation: unclogging of a single pore

dimension = 2;
numPartitions = 64; % Number of partitions in each direction
endTime = 0.8;
dt = 0.5 / numPartitions;
contourPlotLineWidth = 2;
global EPS;

microscaleGrid = FoldedCartesianGrid(dimension, kron(ones(1, dimension), [-0.5, 0.5]), numPartitions*ones(1, dimension));
coord = microscaleGrid.coordinates;
coordCell = mat2cell(coord, ones(1, microscaleGrid.nodes), dimension);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Limiting cylinder
stopRadius = 0.15; % + g.stepSize(1);
stopFun = @(x) stopRadius - abs(x(1));
stopData = cellfun(stopFun, coordCell);
plotStop = microscaleGrid.reshape(stopData);

% lsf = @(x) 0.4 - norm( [2*x(1), x(2)] -[0, -0.5], Inf );
initialLevelSetFunc = stopFun;
initialData = cellfun(initialLevelSetFunc, coordCell);

X = microscaleGrid.reshape(coord(:, 1));
Y = microscaleGrid.reshape(coord(:, 2));
plotInitial = microscaleGrid.reshape(initialData);

exponent = 1;
speedMagnitude = 20;
verticalBorder = 0.1;
borderFun = @(t) verticalBorder + 0.25 * t;
speed = @(t, y) speedMagnitude * ...
    (borderFun(t) - y(2))^exponent * ...
    (borderFun(t) + y(2))^exponent * ...
    (abs(y(2)) < borderFun(t));
% supportBorderInv = 5;
% speed = @(t,y) exp( - ( supportBorderInv * y(2) )^2 / ...
%         ( ( supportBorderInv * y(2) )^2 - 1 ) ) * ...
%         ( abs( y(2) ) < 1/supportBorderInv );

fig = figure;
[~, defHandle] = contour(X, Y, plotInitial, zeros(1, 2), 'LineColor', 'r', ...
    'LineWidth', contourPlotLineWidth);
set(gca, 'xtick', [], 'ytick', []);
axis([-0.5, 0.5, -0.5, 0.5]);
axis equal;
defHandle = copyobj(defHandle, gca);
set(defHandle, 'Zdata', plotStop, 'LineColor', 'g', 'LineWidth', 3);
hand = copyobj(defHandle, gca);
set(hand, 'Zdata', plotInitial, 'LineColor', 'b', 'LineWidth', 3);
%print(fig, '-dpng', 'data/levelSet0.png');
printFileName = 'data/PoreUnclogging/Interfaces_TimeStep';
% print(fig, '-dpng', [ printFileName, num2str(0), '.png' ] );
clear defHandle defHandle2

tic;
% disp( 'Move level set' );
% tic;
%levelSet = solveLevelSetEquation( g, initialData, speed, velocity, ...
%     endTime, dt, 'ContourHandle', hand, 'SaveTime', 'all', 'ContourFigure', fig );
levelSet = solveLevelSetEquationOld(microscaleGrid, initialData, speed, [], ...
    endTime, dt, 'SaveTime', 'all', 'ContourHandle', hand, ...
    'ContourFig', fig, 'Print', false, 'PrintFileName', ...
    printFileName);
% toc
% print(fig, '-dpng', 'data/levelSetEnd.png');

numTimeSlices = size(levelSet, 2);
diffusionTensors = NaN(4, numTimeSlices);
porosities = NaN(numTimeSlices, 1);

%compute diffusion tensor and porosity for each time-step
for ts = 1:numTimeSlices
    disp(['Time step ', num2str(ts)]);

    [A, rhs, isDoF, triangleVolumes] = assembleCellProblem(microscaleGrid, levelSet(:, ts));
    SOL = solveSystemFE(microscaleGrid, A, rhs, isDoF);

    [diffusion, porosity] = computeDiffusionTensor(microscaleGrid, SOL, triangleVolumes);
    diffusionTensors(:, ts) = diffusion(:);
    porosities(ts) = porosity;
end
toc
