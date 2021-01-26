% Test Level-Set equation with modified normal velocity
% --> obatin droplet shape

%parameters
dimension = 2;
numPartitions = 200; % Number of partitions in each direction
dt = 0.5 / numPartitions;
endTime = 0.08;
contourPlotLineWidth = 3;

microscaleGrid = FoldedCartesianGrid(dimension, kron(ones(1, dimension), [-0.5, 0.5]), ...
    numPartitions*ones(1, dimension));
coord = microscaleGrid.coordinates;

levelSetFunc = @(x) 0.35 - norm(x);
coordCell = mat2cell(coord, ones(1, microscaleGrid.nodes), dimension);
initialData = cellfun(levelSetFunc, coordCell);

X = microscaleGrid.reshape(coord(:, 1));
Y = microscaleGrid.reshape(coord(:, 2));
plotInitial = microscaleGrid.reshape(initialData);

% normal interface speed
speed = @(t, x) 0.1; %0 * (  norm( x, Inf ) <= 0.3 );
%variableSpeed = @(t,x) speed(t,x);%* ( ( t <= 0.5*endTime ) - ( t > 0.5*endTime ) );

% artificial exterior velocity field (external flow field)
velocity = [0.1; 0];

toPrint = false;
printFileName = 'data/AnisotropicMicroscale/Test_';
contourFig = figure;
% plotGrid( microscaleGrid, fig );
[~, defHandle] = contour(X, Y, plotInitial, zeros(1, 2), 'LineColor', 'r', ...
    'LineWidth', contourPlotLineWidth);
set(gca, 'DataAspectRatio', [1, 1, 1]);
set(gca, 'xtick', []);
set(gca, 'ytick', []);
hand = copyobj(defHandle, gca);
set(hand, 'Zdata', plotInitial, 'LineColor', 'b', 'LineWidth', contourPlotLineWidth);
if (toPrint)
    print(contourFig, '-depsc', [printFileName, ...
        num2str(0), '.eps']);
end

%advance level-set in time
levelSet = solveLevelSetEquationOld(microscaleGrid, initialData, speed, velocity, ...
    endTime, dt, 'SaveTime', 'all', 'ContourHandle', hand, ...
    'NormalVelocityModifier', @normalVelocityModFunc, 'Reinitialize', false, ...
    'Print', toPrint, 'PrintFileName', printFileName, 'ContourFigure', contourFig);

func = @normalVelocityModFunc;

gridHyPHM = Grid(coord, microscaleGrid.triangles);
[X, Y] = meshgrid(linspace(-0.5, 0.5, numPartitions + 1));
hold on
surf(X, Y, reshape(-double(levelSet(:, end) > 0), numPartitions + 1, numPartitions + 1)', 'EdgeColor', 'None');
colormap gray
% hold on
% gridHyPHM.visualize()


% calculate porosity and diffusivity (should become more unisotropic)
numTimeSlices = size(levelSet, 2);
diffusionTensors = NaN(4, numTimeSlices);
porosities = NaN(numTimeSlices, 1);

for ts = 1:numTimeSlices
    disp(['Time step ', num2str(ts)]);

    [stiffnessMatrix, balanceVector, isDoF, triangleVolumes] = ...
        assembleCellProblem(microscaleGrid, levelSet(:, ts));
    cellProblemSolutions = solveSystemFE(microscaleGrid, stiffnessMatrix, ...
        balanceVector, isDoF);

    [diffusion, porosity] = computeDiffusionTensor(microscaleGrid, ...
        cellProblemSolutions, triangleVolumes);
    diffusionTensors(:, ts) = diffusion(:);
    porosities(ts) = porosity;
end

% normal speed modifier (multiplicative) as a function of the scalar product
% velocity and normal vector
function y = normalVelocityModFunc(x)

positiveLimitValue = 2 / 3;
expFactor = 200;
positive = (x > 0.0);
y = exp(-expFactor*x.^2);
y(positive) = y(positive) + positiveLimitValue ...
    * (1.0 - exp(-expFactor * x(positive).^2));

end
