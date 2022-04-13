% Investigate influence of reinitialization procedure of level-set function

dim = 2;
n = 40; % Number of partitions in each direction
endTime = 0.2;
dt = 0.5 / n;

g = FoldedCartesianGrid(dim, kron(ones(1, dim), [-0.5, 0.5]), n*ones(1, dim));
coord = g.coordinates;

lsf = @(x) 0.2.^40 - norm(x).^40;
coordCell = mat2cell(coord, ones(1, g.nodes), dim);
initialData = cellfun(lsf, coordCell);

X = reshape(g.coordinates(:, 1), g.nodesPerDimension(1), g.nodesPerDimension(2));
Y = reshape(g.coordinates(:, 2), g.nodesPerDimension(1), g.nodesPerDimension(2));
plotInitialData = g.reshape(initialData);

fig = figure;
subplot(2, 3, 1);
surf(X, Y, plotInitialData);
set(gca, 'DataAspectRatio', [1, 1, 1]);
title('Initial level set function');
xlabel('x');
ylabel('y');

subplot(2, 3, 4);
contour(X, Y, plotInitialData, -0.5:0.1:0.5, 'ShowText', 'on');
set(gca, 'DataAspectRatio', [1, 1, 1]);
xlabel('x');
ylabel('y');

speed = @(t, x) -1;
velocity = [];

tic;
phiWithoutReinit = solveLevelSetEquationOld(g, initialData, speed, [], ...
    endTime, dt, 'Reinitialize', false);
toc
plotPhiWithoutReinit = g.reshape(phiWithoutReinit);

tic;
sdf = reinitializeLevelSetOld(g, initialData);
toc
plotSDF = g.reshape(sdf);

subplot(2, 3, 2);
surf(X, Y, plotPhiWithoutReinit);
set(gca, 'DataAspectRatio', [1, 1, 1]);
title('Evolved level set function without reinitialization');
xlabel('x');
ylabel('y');

subplot(2, 3, 5);
contour(X, Y, plotPhiWithoutReinit, -0.5:0.1:0.5, 'ShowText', 'on');
set(gca, 'DataAspectRatio', [1, 1, 1]);
xlabel('x');
ylabel('y');

tic;
phi = solveLevelSetEquationOld(g, initialData, speed, velocity, endTime, dt);
toc
plotPhi = g.reshape(phi);

subplot(2, 3, 3);
surf(X, Y, plotPhi);
set(gca, 'DataAspectRatio', [1, 1, 1]);
title('Evolved level set function with reinitialization');
xlabel('x');
ylabel('y');

subplot(2, 3, 6);
contour(X, Y, plotPhi, -0.5:0.1:0.5, 'ShowText', 'on');
set(gca, 'DataAspectRatio', [1, 1, 1]);
xlabel('x');
ylabel('y');

print(fig, '-dpng', 'reinit.png');