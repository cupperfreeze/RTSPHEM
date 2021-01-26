% Use a few steps of VIIM to generate a non-trivial geometry and calculate
% the smoothed indicator function if \Chi=3

%% Perform evolution of two semincircles via VIIM
% %Full test case
close all
%parameter
StepSizeTime = 0.004;
numTimeSteps = 10;

%setup grid
dimension = 2;
numPartitions = 150;
cellGrid = CartesianGrid(dimension, ...
    kron(ones(1, dimension), [-0.5, 0.5]), ...
    numPartitions*ones(1, dimension));
coord = cellGrid.coordinates;
temp = numPartitions + 1; %temporary variable for plotting
epsilon = 2 / numPartitions;
radius = 0.3;

% Initial indicator function for different subdomains
Xi = ones(size(coord, 1), 1);
Xi((coord(:, 1) > 0) & (sqrt(coord(:, 1).^2 + coord(:, 2).^2) <= radius)) = 2;
Xi((coord(:, 1) <= 0) & (sqrt(coord(:, 1).^2 + coord(:, 2).^2) <= radius)) = 3;

% Velocities at interfaces
% interfaceVelocities(i,j) = velocity between Xi=i and Xi=j
interfaceVelocities = [0, 0, -2; 0, 0, 0; 0, 0, 0];
interfaceVelocities = interfaceVelocities - interfaceVelocities'; %add lower half

%Initial configuration Phi
initialPhiFunc = @(x) max(min(-(radius - epsilon - norm(x, 2)), radius + epsilon - norm(x, 2)), (epsilon - norm(x(1) - 10 * eps, inf))*(abs(x(2)) <= radius)-100*(abs(x(2)) > radius));
coordCell = mat2cell(coord, ones(1, cellGrid.nodes), dimension);
Phi = cellfun(initialPhiFunc, coordCell);

%Calculate distance function d^i only up to this value
restrictDist = 2 * epsilon;
tic
for i = 1:numTimeSteps
    %Calculate distance function from each subdomain
    [distFunctions, Xi] = step1(Xi, Phi, cellGrid, restrictDist);


    %Reconstruct Vonoroi Interface and velocity extension
    [dV, S] = step234(distFunctions, cellGrid, Xi, interfaceVelocities, (10^10 * (mod(i, 10) == 0) + 50)*restrictDist);

    %Update Phi via dV, but only every 10 steps
    if mod(i, 10) == 0
        Phi = epsilon - dV;
    end

    %Evolve Phi
    Phi = levelSetEquationTimeStep(StepSizeTime, 0, Phi, ...
        cellGrid, S', 2);

end
toc


figure
contour(reshape(coord(:, 1), temp, temp), reshape(coord(:, 2), temp, temp), reshape(Phi, temp, temp), linspace(-epsilon, epsilon, 10));
title('Current Phi Contour');

figure
surf(reshape(coord(:, 1), temp, temp), reshape(coord(:, 2), temp, temp), reshape(Phi, temp, temp));
title('Current Phi');

figure
surf(reshape(coord(:, 1), temp, temp), reshape(coord(:, 2), temp, temp), reshape(Xi, temp, temp));
set(gca, 'FontSize', 15) %title('Current Xi');

figure
surf(reshape(coord(:, 1), temp, temp), reshape(coord(:, 2), temp, temp), reshape(dV, temp, temp));
title('Current dV');

figure
surf(reshape(coord(:, 1), temp, temp), reshape(coord(:, 2), temp, temp), reshape(S, temp, temp));
title('Current S');
[interfaceLength, coordTripel] = evaluateInterface(cellGrid, Xi, distFunctions, false)

%% use geometry of last VIIM step and plot smoothed indicator wrt Voronoi extension of solid phases to the whole domain

figure
Chart = Xi;
delta = 0.1;
XiV = Xi;
% for each point in the domain, compute distance to solid phases and
% obtain Voronoi extensions of the sets
for i = 1:numel(Xi)
    x = coord(i, 1);
    y = coord(i, 2);
    dist1 = min(sqrt((x-coord(Xi == 3, 1)).^2 + (y - coord(Xi == 3, 2)).^2));
    dist2 = min(sqrt((x-coord(Xi == 2, 1)).^2 + (y - coord(Xi == 2, 2)).^2));
    XiV(i) = (dist1 > dist2) + 2;
end

% smoothing procedure according to Voronoi extensions
for i = 1:numel(Xi)
    x = coord(i, 1);
    y = coord(i, 2);
    dist1 = min(sqrt((x-coord(XiV ~= 2, 1)).^2 + (y - coord(XiV ~= 2, 2)).^2));
    dist2 = min(sqrt((x-coord(XiV == 2, 1)).^2 + (y - coord(XiV == 2, 2)).^2));
    Chart(i) = (dist1 > delta) * 0 + (dist2 > delta) * 1 + (0 < dist2) * (dist2 <= delta) * (0.5 + dist2 / (2 * delta)) + (0 < dist1) * (dist1 <= delta) * (0.5 - dist1 / (2 * delta));
end

% plot smoothed function
pic = -1 - reshape(Chart, 151, 151);
surf(reshape(coord(:, 1), 151, 151), reshape(coord(:, 2), 151, 151), pic, 'EdgeColor', 'None')
shading interp
set(gca, 'FontSize', 16);
colormap gray
yticks([0, 0.33, 0.4])
yticklabels({'0', '0.33', '0.4'})
ylabel('y');
xlabel('x');
%gca.Color='black';
axis equal
grid off

hold on
[interfaceLength, coordTripel] = evaluateInterface(cellGrid, Xi, distFunctions, false)


% plot smoothed function along y-slices
figure
plot(coord(1:151, 1), pic(:, 75)+2, 'LineWidth', 3);
hold on
plot(coord(1:151, 1), pic(:, 125)+2, 'LineWidth', 3);
hold on
plot(coord(1:151, 1), pic(:, 135)+2, 'LineWidth', 3);
set(gca, 'FontSize', 16);
legend('slice y=0', 'slice y=0.33', 'slice y=0.4')
xlabel('x');
ylabel('(\chi =3)_\delta')
xticks(linspace(-.5, 0.5, 11));
