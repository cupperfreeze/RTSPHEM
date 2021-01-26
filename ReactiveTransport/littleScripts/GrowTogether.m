% Track surface-porosity relation for two circles merging

global EPS
EPS = eps;

speed = -1;
% parameter
stepsize = 1 / 400;
numSteps = 100;
dimension = 2;
numPartitions = 64;

levelSet = nan((numPartitions+1)^2, numSteps);

% setup grid
cellGrid = FoldedCartesianGrid(dimension, ...
    kron(ones(1, dimension), [-0.5, 0.5]), ...
    numPartitions*ones(1, dimension));
coord = cellGrid.coordinates;

initialLevelSetFunc = @(x) max(0.05-norm(x - [0.1, 0], 2), 0.05-norm(x + [0.1, 0], 2));

coordCell = mat2cell(coord, ones(1, cellGrid.nodes), dimension);
levelSet(:, 1) = cellfun(initialLevelSetFunc, coordCell); % initial level-set

% plot initial configuration
[a, b] = meshgrid(-0.5:1/numPartitions:0.5, -0.5:1/numPartitions:0.5);
figure
contourf(a, b, reshape(levelSet(:, 1), numPartitions + 1, numPartitions + 1)', [0, 0], 'linewidth', 2, 'color', 0.5*[1, 1, 1])
axis equal
colormap(gray)
xlabel('x');
ylabel('y');
set(gca, 'FontSize', 15);
clear Surfaces Areas

% Evolve level-set function and track surface area/porosity
for i = 1:numSteps
    [cellProblemSystemMatrix, rhs, isDoF, triangleVolumes, triangleSurfaces] ...
        = assembleCellProblem(cellGrid, levelSet(:, i));
    Surfaces(i, 1) = sum(triangleSurfaces);
    Areas(i, 1) = sum(triangleVolumes);
    levelSet(:, i+1) = levelSetEquationTimeStep(stepsize, 0, levelSet(:, i), cellGrid, speed, 2);
end
[cellProblemSystemMatrix, rhs, isDoF, triangleVolumes, triangleSurfaces] ...
    = assembleCellProblem(cellGrid, levelSet(:, numSteps + 1));
%Surfaces(numSteps+1) =  sum( triangleSurfaces );

% plot terminal configuration
figure
[a, b] = meshgrid(-0.5:1/numPartitions:0.5, -0.5:1/numPartitions:0.5);
contourf(a, b, reshape(levelSet(:, numSteps + 1), numPartitions + 1, numPartitions + 1)', [0, 0], 'linewidth', 2, 'color', 0.5*[1, 1, 1])
axis equal
colormap(gray)
xlabel('x');
ylabel('y');
set(gca, 'FontSize', 15);

% plot graphs of surface/porosity relation and compare to theoretical value
% of two distinct growing circles
figure
x = linspace(0.05, 0.05+numSteps*stepsize, numSteps+1);
Ar = 1 - 2 * pi * x.^2;
Sur = 4 * pi * x;
plot(Ar(Ar > Areas(end)), (Sur(Ar > Areas(end))), 'LineWidth', 2)
hold on
plot((Areas), Surfaces, 'LineWidth', 2)
set(gca, 'FontSize', 15);
legend('reference: disjoint circles', 'merging circles');
xlabel('porosity');
ylabel('surface area');