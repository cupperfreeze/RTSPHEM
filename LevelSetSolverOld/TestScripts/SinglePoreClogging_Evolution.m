%Level-set simulation: clogging of a single pore
% cf. [5] Section 4.1

dimension = 2;
numPartitions = 64; % Number of partitions in each direction
endTime = 1.5; %1 + 10/128;
timeStepSize = 0.5 / numPartitions;

global EPS;

cellGrid = FoldedCartesianGrid(dimension, ...
    kron(ones(1, dimension), [-0.5, 0.5]), ...
    numPartitions*ones(1, dimension));
coord = cellGrid.coordinates;
stopRadius = 0.15; % + g.stepSize(1);

% lsf = @(x) 0.4 - norm( [2*x(1), x(2)] -[0, -0.5], Inf );
initialLevelSetFunc = @(x) max(0.3-norm([2 * x(1), x(2)] - [0, -0.5], 5), ...
    0.3-norm([2 * x(1), x(2)] - [0, 0.5], 5));
coordCell = mat2cell(coord, ones(1, cellGrid.nodes), dimension);
initialData = cellfun(initialLevelSetFunc, coordCell);

X = cellGrid.reshape(coord(:, 1));
Y = cellGrid.reshape(coord(:, 2));
clear coord;

plotInitial = cellGrid.reshape(initialData);

stopFun = @(x) stopRadius - abs(x(1));
stopData = cellfun(stopFun, coordCell);
clear coordCell;

plotStop = cellGrid.reshape(stopData);

speed = @(t, x) -0.2 * (abs(x(1)) < stopRadius) + 0.2 * (abs(x(1)) >= stopRadius);

%velocity = zeros( cellGrid.nodes, dimension );
%velocity=[1,0]'*0.8;
velocity = [];

fig = figure;
[~, defHandle] = contour(X, Y, plotInitial, zeros(1, 2), 'LineColor', 'r', 'LineWidth', 3);
% set( gca, 'DataAspectRatio', [1 1 1]);
% set( gca, 'fontsize', 20 );
set(gca, 'xtick', [], 'ytick', []);
axis([-0.5, 0.5, -0.5, 0.5])
axis equal
% xlabel('x');
% ylabel('y');
% handle = copyobj( defHandle, gca );
defHandle2 = copyobj(defHandle, gca);
set(defHandle2, 'Zdata', plotStop, 'LineColor', 'g', 'LineWidth', 3);
hand = copyobj(defHandle2, gca);
set(hand, 'Zdata', plotInitial, 'LineColor', 'b', 'LineWidth', 3);
%print(fig, '-dpng', 'data/levelSet0.png');
dataFolder = 'data/SinglePoreClogging/';
printFileName = [dataFolder, 'TimeStep_'];
toPrint = false;
% print(fig, '-dpng', [ printFileName, num2str(0), '.png' ] );
if (toPrint)
    print(fig, '-depsc', [printFileName, num2str(0), '.eps']);
end

tic;
% disp( 'Move level set' );
% tic;
%levelSet = solveLevelSetEquation( g, initialData, speed, velocity, ...
%     endTime, dt, 'ContourHandle', hand, 'SaveTime', 'all', 'ContourFigure', fig );
levelSet = solveLevelSetEquationOld(cellGrid, initialData, speed, velocity, ...
    endTime, timeStepSize, 'SaveTime', 'all', 'Reinitialize', false, 'ContourHandle', hand, ...
    'ContourFig', fig, 'Print', toPrint, 'PrintFileName', ...
    printFileName);
% toc
% print(fig, '-dpng', 'data/levelSetEnd.png');

clear defHandle defHandle2 fig hand;

numTimeSlices = size(levelSet, 2);
diffusionTensors = NaN(4, numTimeSlices);
porosities = NaN(numTimeSlices, 1);

%compute diffusion tensor and porosity for each time-step
for ts = 1:numTimeSlices
    disp(['Time step ', num2str(ts)]);

    [A, rhs, isDoF, triangleVolumes] = assembleCellProblem(cellGrid, levelSet(:, ts));
    SOL = solveSystemFE(cellGrid, A, rhs, isDoF);

    [diffusion, porosity] = computeDiffusionTensor(cellGrid, SOL, triangleVolumes);
    diffusionTensors(:, ts) = diffusion(:);
    porosities(ts) = porosity;
end
toc

clear diffusion porosity;

%dateString = datestr( now, 'yyyymmdd_HHMM' );
%save( [ dataFolder, 'Data_', dateString, '.mat' ] );
