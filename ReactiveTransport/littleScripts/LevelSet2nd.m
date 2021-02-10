% Script to evaluate the accuracy of 1st and 2nd order level-set solver

global EPS
EPS = eps;
figure
steps = 50;
outstep = 50;
radius = 0.3;
speed = 1;

numPartitions = 64;

dimension = 2;
stepsize = 0.005;
cellGrid = FoldedCartesianGrid(dimension, ...
    kron(ones(1, dimension), [-0.5, 0.5]), ...
    numPartitions*ones(1, dimension));
LEVELS1 = cell(steps, 1);
LEVELS2 = cell(steps, 1);

coord = cellGrid.coordinates;


initialLevelSetFunc = @(x) radius - norm([x(1), x(2)]-[0.0, 0.0], 2); %radius-norm(x,2)-norm(x-[0.2,0.1],3);
coordCell = mat2cell(coord, ones(1, cellGrid.nodes), dimension);
levelSet = cellfun(initialLevelSetFunc, coordCell);
LEVELS1{1} = levelSet;
LEVELS2{1} = levelSet;
tic
for i = 1:(steps - 1)
    LEVELS1{i+1} = levelSetEquationTimeStep(stepsize, 0, LEVELS1{i}, cellGrid, speed, 1);
end
toc

tic
for i = 1:(steps - 1)
    LEVELS2{i+1} = levelSetEquationTimeStep(stepsize, 0, LEVELS2{i}, cellGrid, speed);
end
toc
l = -0.5:1 / numPartitions:0.5;
[a, b] = meshgrid(l, l);

contour(a, b, reshape(LEVELS1{1}, numPartitions + 1, numPartitions + 1), [0, 0], 'linewidth', 2, 'color', 'y')
hold on
contour(a, b, reshape(LEVELS1{outstep}, numPartitions + 1, numPartitions + 1), [0, 0], 'linewidth', 3, 'color', 'r');
hold on
contour(a, b, reshape(LEVELS2{outstep}, numPartitions + 1, numPartitions + 1), [0, 0], 'linewidth', 2, 'color', 'b');
hold on
initialLevelSetFunc = @(x) radius - speed * stepsize * (outstep - 1) - norm([x(1), x(2)]-[0.0, 0.0], 2); %radius-norm(x,2)-norm(x-[0.2,0.1],3);
coordCell = mat2cell(coord, ones(1, cellGrid.nodes), dimension);
levelSet = cellfun(initialLevelSetFunc, coordCell);
%contour(a,b,reshape(levelSet,numPartitions+1, numPartitions+1),[0,0], 'linewidth', 2,'color','k');
set(gca, 'FontSize', 15)
axis equal
xlabel('x')
ylabel('y')
legend('initial', '1st order', '2nd order', 'analytical')
%title('Level-set solver accuracy');
