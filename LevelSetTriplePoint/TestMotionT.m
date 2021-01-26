%Full test case for initial T-shape
close all
%parameter
StepSizeTime = 0.0025;
numTimeSteps = 40;

%setup grid
dimension = 2;
numPartitions = 100;
cellGrid = CartesianGrid(dimension, ...
    kron(ones(1, dimension), [-0.5, 0.5]), ...
    numPartitions*ones(1, dimension));
coord = cellGrid.coordinates;
temp = numPartitions + 1; %temporary variable for plotting

% Initial indicator function for different subdomains
Xi = ones(size(coord, 1), 1);
Xi((coord(:, 2) < 0) & (coord(:, 1) <= 0)) = 2;
Xi((coord(:, 2) < 0) & (coord(:, 1) > 0)) = 3;

% Velocities at interfaces
% interfaceVelocities(i,j) = velocity between Xi=i and Xi=j
interfaceVelocities = [0, 0, -2; 0, 0, 0; 0, 0, 0];
interfaceVelocities = interfaceVelocities - interfaceVelocities'; %add lower half

%Initial configuration Phi
epsilon = 2 / numPartitions;
initialPhiFunc = @(x) max((epsilon-norm([x(1) + eps], 2))*(x(2) < 0)-10*(x(2) >= 0), epsilon-norm([x(2) + eps], 2));
coordCell = mat2cell(coord, ones(1, cellGrid.nodes), dimension);
Phi = cellfun(initialPhiFunc, coordCell);

%Calculate distance function d^i only up to this value
restrictDist = 2 * epsilon;
tic
for i = 1:numTimeSteps
    %Calculate distance function from each subdomain
    [distFunctions, Xi] = step1(Xi, Phi, cellGrid, restrictDist);


    %Reconstruct Vonoroi Interface and velocity extension
    [dV, S] = step234(distFunctions, cellGrid, Xi, interfaceVelocities, 0.1);

    %Update Phi, but only every 10 steps
    if mod(i, 10) == 0
        Phi = epsilon - dV;
    end

    %Evolve Phi
    Phi = levelSetEquationTimeStep(StepSizeTime, 0, Phi, ...
        cellGrid, S', 1);

end
toc
%plot evolved functions
[distFunctions, Xi] = step1(Xi, Phi, cellGrid, restrictDist);
[interfaceLength, coordTripel] = evaluateInterface(cellGrid, Xi, distFunctions, false)


figure
contour(reshape(coord(:, 1), temp, temp), reshape(coord(:, 2), temp, temp), reshape(Phi, temp, temp), linspace(-epsilon, epsilon, 10));
title('Current Phi Contour');

figure
surf(reshape(coord(:, 1), temp, temp), reshape(coord(:, 2), temp, temp), reshape(Phi, temp, temp));
title('Current Phi');

figure
surf(reshape(coord(:, 1), temp, temp), reshape(coord(:, 2), temp, temp), reshape(Xi, temp, temp));
title('Current Xi');

figure
surf(reshape(coord(:, 1), temp, temp), reshape(coord(:, 2), temp, temp), reshape(dV, temp, temp));
title('Current dV');

figure
surf(reshape(coord(:, 1), temp, temp), reshape(coord(:, 2), temp, temp), reshape(S, temp, temp));
title('Current S');