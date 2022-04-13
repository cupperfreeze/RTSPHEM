% Level - set simulation, non-uniform expansion and contraction of a circle

dim = 2;
n = 80; % Number of partitions in each direction
endTime = 0.6;
dt = 0.5 / n;

g = FoldedCartesianGrid(dim, kron(ones(1, dim), [-0.5, 0.5]), n*ones(1, dim));
coord = g.coordinates;

lsf = @(x) 0.2 - norm(x-[0.25, 0]);
coordCell = mat2cell(coord, ones(1, g.nodes), dim);
initialData = cellfun(lsf, coordCell);

X = g.reshape(coord(:, 1));
Y = g.reshape(coord(:, 2));
plotInitial = g.reshape(initialData);

speed = @(t, x) -35 * (x(1) - 0.3).^2 * (x(1) + 0.5).^2 * (x(1) <= 0.3);
variableSpeed = @(t, x) speed(t, x) * ((t <= 0.5 * endTime) - (t > 0.5 * endTime));
velocity = [];

fig = figure;
[~, defHandle] = contour(X, Y, plotInitial, zeros(1, 2), 'LineColor', 'r');
set(gca, 'DataAspectRatio', [1, 1, 1]);
xlabel('x');
ylabel('y');
handle = copyobj(defHandle, gca);
set(handle, 'Zdata', plotInitial, 'LineColor', 'b');
%print(fig, '-dpng', 'levelSet0.png');

phi = solveLevelSetEquationOld(g, initialData, variableSpeed, ...
    velocity, endTime, dt, 'ContourHandle', handle, 'SaveTime', 'all', ...
    'ContourFigure', fig);