% Test of XFEM for auxiliary diffusion cell problems.
% Computes cell problem solutions for a fixed interface defined implicitly
% via a level set function and plots them.
contourPlotLineWidth = 1;
global EPS;

numPartitions = 200;
microscaleGrid = FoldedCartesianGrid(2, [-0.5, 0.5, -0.5, 0.5], numPartitions*ones(1, 2));
initialLevelSetFunc = @(x, y) 0.21 - norm([x, 0.8 * y], 2);
levelSet = arrayfun(initialLevelSetFunc, ...
    microscaleGrid.coordinates(:, 1), microscaleGrid.coordinates(:, 2));
% lsf = @(x) max( ...
%     min( min( x(1) + 0.25, x(2) + 0.25 ), 0.25 - norm( x + 0.25, 1 ) ), ...
%     min( min( -(x(1) - 0.25), -(x(2) - 0.25) ), 0.25 - norm( -(x - 0.25), 1 ) ) );
% coord = g.coordinates;
% coordCell = mat2cell( coord, ones(1,g.nodes), 2 );
% levelSet = cellfun( lsf, coordCell );

disp('Calculating cell problem solutions...');
tic;

[A, rhs, isDoF, triangleVolumes] = assembleCellProblem(microscaleGrid, levelSet);
cellProblemSolution = solveSystemFE(microscaleGrid, A, rhs, isDoF);
[diffusion, porosity] = computeDiffusionTensor(microscaleGrid, ...
    cellProblemSolution, triangleVolumes);

simulationTime = toc;
disp(['    ...done in ', num2str(simulationTime), ' seconds.']);

X = microscaleGrid.reshape(microscaleGrid.coordinates(:, 1));
Y = microscaleGrid.reshape(microscaleGrid.coordinates(:, 2));

figure;
% Plot solution of first cell problem solution together with the interface.
subplot(2, 2, 1);
h = trisurf(microscaleGrid.triangles, microscaleGrid.coordinates(:, 1), ...
    microscaleGrid.coordinates(:, 2), cellProblemSolution(:, 1),'EdgeColor','none');
axis equal;
xlabel('x');
ylabel('y');
hold on
% Assert that the contour is plotted above all values of the solution so
% that it is visible completely in 2D view.
elevation(1) = max(cellProblemSolution(:, 1)) + EPS;
contour3(X, Y, microscaleGrid.reshape(levelSet + elevation(1)), ...
    elevation(1)*ones(2, 1), 'LineColor', 'k', 'LineWidth', contourPlotLineWidth);
view(2); % 2D view

% Plot solution of second cell problem solution together with the interface.
subplot(2, 2, 2);
trisurf(microscaleGrid.triangles, microscaleGrid.coordinates(:, 1), ...
    microscaleGrid.coordinates(:, 2), cellProblemSolution(:, 2),'EdgeColor','none');
axis equal;
xlabel('x');
ylabel('y');
hold on
% Assert that the contour is plotted above all values of the solution so
% that it is visible completely in 2D view.
elevation(2) = max(cellProblemSolution(:, 2)) + eps;
contour3(X, Y, microscaleGrid.reshape(levelSet + elevation(2)), ...
    elevation(2)*ones(2, 1), 'LineColor', 'k', 'LineWidth', contourPlotLineWidth);
view(2); % 2D view

% 3D plot of the cell problem solutions.
subplot(2, 2, 3);
trisurf(microscaleGrid.triangles, microscaleGrid.coordinates(:, 1), ...
    microscaleGrid.coordinates(:, 2), cellProblemSolution(:, 1),'EdgeColor','none');
axis equal;
xlabel('x');
ylabel('y');
subplot(2, 2, 4);
trisurf(microscaleGrid.triangles, microscaleGrid.coordinates(:, 1), ...
    microscaleGrid.coordinates(:, 2), cellProblemSolution(:, 2),'EdgeColor','none');
axis equal;
xlabel('x');
ylabel('y');
