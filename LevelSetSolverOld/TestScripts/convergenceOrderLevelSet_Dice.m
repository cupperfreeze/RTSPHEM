% Test convergence order of level-set solver in case of rectangles

dim = 2;
initialExponent = 0; % Number of partitions in each direction
ratio = 2;
n = ratio^initialExponent;
numIterations = 7;
dt = 0.5 / (8 * ratio^(initialExponent + numIterations));
endTime = 0.1;

disp(['h = ', num2str(1 / n)]);

grids = cell(1+numIterations, 1);
grids{1} = FoldedCartesianGrid(dim, kron(ones(1, dim), [-0.5, 0.5]), n*ones(1, dim));
coord = grids{1}.coordinates;

lsf = @(x) max( ...
    max( ...
    max( ...
    max( ...
    0.1 - norm(x-[0.25, 0.25]), ...
    0.1 - norm(x-[-0.25, -0.25])), ...
    0.1 - norm(x)), ...
    0.1 - norm(x - [0.25, -0.25])), ...
    0.1-norm(x - [-0.25, 0.25]));
coordCell = mat2cell(coord, ones(1, grids{1}.nodes), dim);
initialData = cellfun(lsf, coordCell);

Xcoarse = grids{1}.reshape(coord(:, 1));
Ycoarse = grids{1}.reshape(coord(:, 2));
plotInitial = grids{1}.reshape(initialData);

speed = @(t, x) -0.5;
%speed = @(t,x) -1 * ( abs( x(1) + 0.25 ) < 0.1 | abs( x(1) - 0.25 ) < 0.1 );
%variableSpeed = @(t,x) speed(t,x);%* ( ( t <= 0.5*endTime ) - ( t > 0.5*endTime ) );

velocity = [];

% fig = figure;
% [~, defHandle] = contour( Xcoarse, Ycoarse, plotInitial, zeros(1,2), 'LineColor', 'r' );
% set(gca,'DataAspectRatio',[1 1 1]);
% xlabel('x');
% ylabel('y');
% handle = copyobj( defHandle, gca );
% set( handle, 'Zdata', plotInitial, 'LineColor', 'b' );
%print(fig, '-dpng', 'levelSet0.png');

levelSet = cell(1+numIterations, 1);

levelSet{1} = solveLevelSetEquationOld(grids{1}, initialData, speed, velocity, ...
    endTime, dt, 'SaveTime', 'last'); %, 'ContourFig', fig, 'ContourHandle', handle );

for i = 1:numIterations
    n = ratio * n;

    disp(['h = ', num2str(1 / n)]);

    disp('Preprocessing...');
    tic;

    grids{1+i} = FoldedCartesianGrid(dim, kron(ones(1, dim), [-0.5, 0.5]), ...
        n*ones(1, dim));
    coord = grids{1+i}.coordinates;

    coordCell = mat2cell(coord, ones(1, grids{1 + i}.nodes), dim);
    initialData = cellfun(lsf, coordCell);

    velocity = [];

    %     figFine = figure;
    %     [~, defHandle] = contour( Xfine, Yfine, plotInitial, zeros(1,2), 'LineColor', 'r' );
    %     set(gca,'DataAspectRatio',[1 1 1]);
    %     xlabel('x');
    %     ylabel('y');
    %     handleFine = copyobj( defHandle, gca );
    %     set( handleFine, 'Zdata', plotInitial, 'LineColor', 'b' );
    %print(fig, '-dpng', 'levelSet0.png');

    toc
    disp('Level set evolution...');
    tic;

    levelSet{1+i} = solveLevelSetEquationOld(grids{1 + i}, initialData, speed, velocity, ...
        endTime, dt, 'SaveTime', 'last');

    toc

end

l1Error = zeros(numIterations, 1);
l2Error = zeros(numIterations, 1);
lInfError = zeros(numIterations, 1);
for i = 1:numIterations
    l1Error(i) = L1Error(grids{i}, levelSet{i}, grids{i + 1}, levelSet{i + 1});
    l2Error(i) = L2Error(grids{i}, levelSet{i}, grids{i + 1}, levelSet{i + 1});
    lInfError(i) = LInfError(grids{i}, levelSet{i}, grids{i + 1}, levelSet{i + 1});
end

EOC_l1 = zeros(numIterations-1, 1);
EOC_l2 = zeros(numIterations-1, 1);
EOC_lInf = zeros(numIterations-1, 1);
for j = 1:numIterations - 1
    EOC_l1(j) = log2(l1Error(j)/l1Error(j + 1));
    EOC_l2(j) = log2(l2Error(j)/l2Error(j + 1));
    EOC_lInf(j) = log2(lInfError(j)/lInfError(j + 1));
end