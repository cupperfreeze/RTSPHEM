dim = 2;
n = 64; % Number of partitions in each direction
endTime = 0.2;
dt = 0.5 / n;

g = FoldedCartesianGrid(dim, kron(ones(1, dim), [-0.5, 0.5]), n*ones(1, dim));
coord = g.coordinates;

lsf = @(x) max(max(max(max(0.1 - norm(x-[-0.25, -0.25], Inf), 0.1 - norm(x-[0.25, 0], Inf)), ...
    0.1 - norm(x - [-0.25, 0.25], Inf)), 0.1 - norm(x - [0.25, 0.5], Inf)), ...
    0.1-norm(x - [0.25, -0.5], Inf));
coordCell = mat2cell(coord, ones(1, g.nodes), dim);
initialData = cellfun(lsf, coordCell);

X = g.reshape(coord(:, 1));
Y = g.reshape(coord(:, 2));
plotInitial = g.reshape(initialData);

speed = @(t, x) 0.25;
%speed = @(t,x) -1 * ( abs( x(1) + 0.25 ) < 0.1 | abs( x(1) - 0.25 ) < 0.1 );
%variableSpeed = @(t,x) speed(t,x);%* ( ( t <= 0.5*endTime ) - ( t > 0.5*endTime ) );

velocity = [];

fig = figure;
[~, defHandle] = contour(X, Y, plotInitial, zeros(1, 2), 'LineColor', 'r');
set(gca, 'DataAspectRatio', [1, 1, 1]);
xlabel('x');
ylabel('y');
handle = copyobj(defHandle, gca);
set(handle, 'Zdata', plotInitial, 'LineColor', 'b');
%print(fig, '-dpng', 'levelSet0.png');

phi = solveLevelSetEquationOld(g, initialData, speed, velocity, ...
    endTime, dt, 'ContourHandle', handle, 'ContourFig', fig);
