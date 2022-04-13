% Test of microscale coupling between level set evolution and calculation
% of effective parameters with triangle interfaces.
% Simulates shrinking of triangles and the resulting change in effective
% parameters.
EPS = eps
contourPlotLineWidth = 1;
dimension = 2;
numPartitions = 64; % Number of partitions in each direction
dt = 0.5 / numPartitions;
numTimeSteps = 10;
endTime = 0.3;

microscaleGrid = FoldedCartesianGrid(dimension, kron(ones(1, dimension), [-0.5, 0.5]), numPartitions*ones(1, dimension));
coord = microscaleGrid.coordinates;

%lsf = @(x) max( ...
%    min( min( x(1) + 1/3, x(2) + 1/3 ), 1 - norm( x + 1/3, 1 ) ), ...
%    min( min( -(x(1) - 1/3), -(x(2) - 1/3) ), 1/3 - norm( -(x - 1/3), 1 ) ) );
triangleEdgeLength = sqrt(1-0.2);
triangleFunc = @(x) max(max(-x(1) - (1 / 6 + triangleEdgeLength / 3), ...
    -x(2) - (1 / 6 + triangleEdgeLength / 3)), ...
    (x(1) + x(2))-(triangleEdgeLength / 3 - 1 / 3));
initialLevelSetFunc = @(x) -min(triangleFunc(x), triangleFunc(-x));
coordCell = mat2cell(coord, ones(1, microscaleGrid.nodes), dimension);
initialData = cellfun(initialLevelSetFunc, coordCell);

X = microscaleGrid.reshape(coord(:, 1));
Y = microscaleGrid.reshape(coord(:, 2));
plotInitial = microscaleGrid.reshape(initialData);

speed = @(t, x) (sqrt(0.8) - sqrt(0.2));

fig = figure;
%plotGrid( microscaleGrid, fig );
[~, defHandle] = contour(X, Y, plotInitial, zeros(1, 2), 'LineColor', 'r', ...
    'LineWidth', contourPlotLineWidth);
set(gca, 'DataAspectRatio', [1, 1, 1]);
xlabel('x');
ylabel('y');
hand = copyobj(defHandle, gca);
set(hand, 'Zdata', plotInitial, 'LineColor', 'b');

tic;
t_ = 0;
phi = initialData(:);
phiOld = phi;
levelSet = NaN(numel(initialData), ceil(endTime / dt)+1);
levelSet(:, 1) = phi;

timeStep = 1;

diffusionTensors = cell(1, ceil(endTime / dt)+1);
porosities = cell(ceil(endTime / dt)+1, 1);

[A, rhs, isDoF, triangleVolumes] = assembleCellProblem(microscaleGrid, levelSet(:, timeStep));
SOL = solveSystemFE(microscaleGrid, A, rhs, isDoF);
[diffusion, porosity] = computeDiffusionTensor(microscaleGrid, SOL, triangleVolumes);
diffusionTensors{timeStep}(:, 1) = diffusion(:);
porosities{timeStep}(1) = porosity;
while (t_ <= endTime + EPS)

    assert(all(phi >= phiOld - EPS));

    if (t_ + dt > endTime - EPS)
        dt = endTime - t_;
        if (abs(dt) < EPS)
            break;
        end
    end
    timeCell = num2cell(t_*ones(microscaleGrid.nodes, 1));
    normalSpeed = cellfun(speed, timeCell, coordCell);

    phi = levelSetEquationTimeStep(t_+dt, t_, phi, microscaleGrid, ...
        normalSpeed);

    plotLS(microscaleGrid, phi, timeStep, hand, fig);
    t_ = t_ + dt;
    phiOld = phi;
    timeStep = timeStep + 1;
    levelSet(:, timeStep) = phi;

    [A, rhs, isDoF, triangleVolumes] = assembleCellProblem(microscaleGrid, levelSet(:, timeStep));
    SOL = solveSystemFE(microscaleGrid, A, rhs, isDoF);
    [diffusion, porosity] = computeDiffusionTensor(microscaleGrid, SOL, triangleVolumes);
    diffusionTensors{timeStep}(:, 1) = diffusion(:);
    porosities{timeStep}(1) = porosity;

end % while
levelSet(:, all(isnan(levelSet))) = [];
toc


function plotLS(grid, phi, timeStep, handle, fig)
contourPlotLineWidth = 1;
dim = grid.dimension;
if (dim < 2)
    PHI = phi;
else
    PHI = grid.reshape(phi);
end

set(handle, 'Zdata', PHI, 'LineColor', 'b', 'LineWidth', contourPlotLineWidth);
%     if ( ~isempty( fig ) )
%         print( fig, '-dpng', [ 'data/levelSet', num2str( timeStep ), '.png' ] );
%     end
drawnow;
end
