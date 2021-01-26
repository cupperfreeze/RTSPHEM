% Evaluate order of convergence for Voronoi Implicit Interface Method
% for normal advection
% results of this script are saved in Data file 'LevelConv2eps.mat'

k = 7; %number of refinement levels
containerXi = cell(k, 1);
containerPhi = cell(k, 1);
containerDistFunctions = cell(k, 1);
containerInterfaceLength = cell(k, 1);
containerCoordTripel = cell(k, 1);
containerTime = cell(k, 1);

for l = 1:k
    %Full test case
    close all
    %parameter
    StepSizeTime = 1 / 40 / 2^(l - 1);
    numTimeSteps = 10 * 2^(l - 1);

    %setup grid
    dimension = 2;
    numPartitions = 20 * 2^(l - 1);
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
    interfaceVelocities = [0, -0, -1; 0, 0, 0; 0, 0, 0];
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
        [dV, S] = step234(distFunctions, cellGrid, Xi, interfaceVelocities, 15/numPartitions);

        %Update Phi, but only every 10 steps
        if mod(i, 10) == 0
            Phi = epsilon - dV;
        end

        %Evolve Phi
        Phi = levelSetEquationTimeStep(StepSizeTime, 0, Phi, ...
            cellGrid, S', 2);

    end
    containerTime{l} = toc;
    containerTime{l}
    %plot evolved functions
    [distFunctions, Xi] = step1(Xi, Phi, cellGrid, restrictDist);
    [interfaceLength, coordTripel] = evaluateInterface(cellGrid, Xi, distFunctions, false)


    containerXi{l} = Xi;
    containerPhi{l} = Phi;
    containerDistFunctions{l} = distFunctions;
    containerInterfaceLength{l} = interfaceLength;
    containerCoordTripel{l} = coordTripel;
    hold on
    Circle = 0.25 - sqrt(coord(:, 1).^2+coord(:, 2).^2);
    contour(reshape(coord(:, 1), temp, temp), reshape(coord(:, 2), temp, temp), reshape(Circle, temp, temp), [-0.1, 0, 1]);
end

figure
plot(1:k, log2([containerTime{:}]), 'LineWidth', 3)
hold on
plot(1:k, 2*(1:k)-1, 'LineWidth', 3)
legend('data', 'quadratic rate')
title('time');
xlabel('refinement level');
ylabel('log2 of time in s');
set(gca, 'FontSize', 15);

figure
lul = [containerCoordTripel{:}];
plot(1:k, log2(abs(-0.25 - lul(1:2:2 * k))), 'LineWidth', 3);
%plot(1:k,log2(sqrt((-sqrt(0.75)*0.25-lul(1:2:2*k)).^2)+(lul(2:2:2*k)-0.125).^2),'LineWidth',3);
hold on
plot(1:k, -(1:k)-1, 'LineWidth', 3)
legend('data', 'linear rate')
title('Triple point convergence');
xlabel('refinement level');
ylabel('log2 distance');
set(gca, 'FontSize', 15);