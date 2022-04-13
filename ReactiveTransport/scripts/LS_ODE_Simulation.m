%%Reactive ODE case, level-set simulation
% cf. [0], Section 4.3


close all
clear all
%parameters
tic
currentTime = 0;
numTimeSteps = 0;
finishTime = 2;
saveTime = 0;

%setup grid
dimension = 2;
numPartitions = 200;
cellGrid = CartesianGrid(dimension, ...
    kron(ones(1, dimension), [-0.5, 0.5]), ...
    numPartitions*ones(1, dimension));
coord = cellGrid.coordinates;
temp = numPartitions + 1; %temporary variable for plotting
epsilon = 2 / numPartitions;
levelSave = zeros((numPartitions+1)^2, 200, 'single');

% Initial indicator function for different subdomains
Xi = ones(size(coord, 1), 1);
Xi((coord(:, 1) >= 0) & (sqrt(coord(:, 1).^2 + coord(:, 2).^2) <= 0.2)) = 2;
Xi((coord(:, 1) < 0) & (sqrt(coord(:, 1).^2 + coord(:, 2).^2) <= 0.2)) = 3;


% Velocities at interfaces
% interfaceVelocities(i,j) = velocity between Xi=i and Xi=j
InitialConcentration = [2; 1; 1];
saveC = InitialConcentration';
interfaceVelocities = [0, 0, 0; 0, 0, 0; 0, 0, 0];
interfaceVelocities = interfaceVelocities - interfaceVelocities'; %add lower half

%Initial configuration Phi
initialPhiFunc = @(x) max(min(-(0.2 - epsilon - norm(x, 2)), 0.2 + epsilon - norm(x, 2)), (epsilon - norm(x(1), inf))*(abs(x(2)) <= 0.2)-100*(abs(x(2)) > 0.2));
coordCell = mat2cell(coord, ones(1, cellGrid.nodes), dimension);
Phi = cellfun(initialPhiFunc, coordCell);

%Calculate distance function d^i only up to this value
restrictDist = 3 * epsilon;

%Calculate distance function from each subdomain
distFunctions = step1(Xi, Phi, cellGrid, restrictDist);


interfaceLength{numTimeSteps+1} = [0.2 * pi, 0.2 * pi, 0.4];
coordTripel{numTimeSteps+1} = [0, -0.2; 0, 0.2];

c = InitialConcentration;
saveSolidArea = 0.2^2 * pi;
saveDArea = saveSolidArea / 2;
[distFunctions, Xi] = step1(Xi, Phi, cellGrid, restrictDist);

while currentTime < finishTime - eps
    %Calculate distance function from each subdomain

    numTimeSteps = numTimeSteps + 1;
    interfaceVelocities = [0, -(c(2) / c(1) - 1), -(c(2) * c(3) - 1); 0, 0, 0; 0, 0, 0];
    interfaceVelocities = interfaceVelocities - interfaceVelocities';

    CFL = min(0.1/max(interfaceVelocities(:))/numPartitions, inf);
    StepSizeTime = min(CFL, finishTime-currentTime);
    currentTime = currentTime + StepSizeTime;
    saveTime = [saveTime; currentTime];


    %Reconstruct Vonoroi Interface and velocity extension
    [dV, S] = step234(distFunctions, cellGrid, Xi, interfaceVelocities, 0.25);

    %Update Phi via dV, but only every 10 steps
    if mod(numTimeSteps, 40) == 0
        Phi = epsilon - dV;
    end

    %Evolve Phi

    Phi = levelSetEquationTimeStep(StepSizeTime, 0, Phi, ...
        cellGrid, S', 1);

    [distFunctions, Xi] = step1(Xi, Phi, cellGrid, restrictDist);
    [dV, S] = step234(distFunctions, cellGrid, Xi, interfaceVelocities, 0.15);
    levelSet = dV;
    levelSet(Xi == 1) = -levelSet(Xi == 1);
    levelSet(Xi == 2) = max(levelSet(Xi == 2), 10^(-4));
    levelSet(Xi == 3) = max(levelSet(Xi == 3), 10^(-4));

    levelSave(:, numTimeSteps) = levelSet;

    %compute mineral volumes and interface lengths
    [~, ~, ~, triangleVolumes, ~] ...
        = assembleCellProblem(cellGrid, levelSet);
    saveSolidArea = [saveSolidArea; 1 - sum(triangleVolumes)];

    levelTemp = dV;
    levelTemp(Xi ~= 2) = -levelTemp(Xi ~= 2);
    [~, ~, ~, triangleVolumes, ~] ...
        = assembleCellProblem(cellGrid, levelTemp);
    saveDArea = [saveDArea; 1 - sum(triangleVolumes)];

    [interfaceLength{numTimeSteps + 1}, coordTripel{numTimeSteps + 1}] = evaluateInterface(cellGrid, Xi, distFunctions, true); %, numTimeSteps, 1/StepSizeTime/10);

    % ODE system for solute concentrations, discretized by implicit Euler
    f = @(x) (x - c) / StepSizeTime - 1 / (1 - saveSolidArea(end)) * ...
        [interfaceLength{numTimeSteps + 1}(1) * (x(2) / x(1) - 1) * (20 + x(1)) + interfaceLength{numTimeSteps + 1}(2) * (x(2) * x(3) - 1) * x(1); ...
        -interfaceLength{numTimeSteps + 1}(2) * (x(2) * x(3) - 1) * (4 - x(2)) - interfaceLength{numTimeSteps + 1}(1) * (x(2) / x(1) - 1) * (20 - x(2)); ...
        -interfaceLength{numTimeSteps + 1}(2) * (x(2) * x(3) - 1) * (4 - x(3)) + interfaceLength{numTimeSteps + 1}(1) * (x(2) / x(1) - 1) * x(3)];

    % solve ODE
    options = optimoptions('fsolve', 'Display', 'none', 'OptimalityTolerance', 10^(-7));
    [c, ~, flag] = fsolve(f, c, options);
    assert(abs(1 - flag) < eps, 'Nonlinear solver failed');
    saveC = [saveC; c'];
end

%% Make Plots

%plot volumes over time
subplot(3, 1, 1)
plot(saveTime, saveDArea, 'LineWidth', 3)
hold on
plot(saveTime, saveSolidArea-saveDArea, 'LineWidth', 3)
hold on
plot(saveTime, repmat(0.0339, numel(saveTime), 1));
hold on
plot(saveTime, repmat(0.1333, numel(saveTime), 1));
set(gca, 'fontsize', 15);
legend('volume D', 'volume P', 'Location', 'northwest');
%title('solid area');
ylabel('Volume')
xlabel('t');
ylim([0, 0.16])

conservation(:, 1) = saveC(:, 1) .* (1 - saveSolidArea) - saveDArea * 20;
conservation(:, 2) = saveC(:, 2) .* (1 - saveSolidArea) + saveDArea * 20 + (saveSolidArea - saveDArea) * 4;
conservation(:, 3) = saveC(:, 3) .* (1 - saveSolidArea) + (saveSolidArea - saveDArea) * 4;
conservation = conservation ./ conservation(1, :);

% plot mass conservation over time
subplot(3, 1, 2)
plot(saveTime, conservation(:, 1), 'LineWidth', 3)
hold on
plot(saveTime, conservation(:, 2), 'LineWidth', 3)
hold on
plot(saveTime, conservation(:, 3), 'LineWidth', 3)
legend('A', 'B', 'C')
set(gca, 'fontsize', 15);
xlabel('t');
ylabel('M_{rel}');
ylim([0.85, 1.1])


%plot surface area over time
subplot(3, 1, 3)
temp2 = cell2mat(interfaceLength);
plot(saveTime, temp2(1:3:end), 'LineWidth', 3)
hold on
plot(saveTime, temp2(2:3:end), 'LineWidth', 3)
hold on
legend('surface D', 'surface P', 'Location', 'northwest');
set(gca, 'fontsize', 15);
xlabel('t');
ylabel('Surface')
ylim([0, 1.5])


% Plot concentrations over time
figure
plot(saveTime, saveC(:, 1), '-*', 'MarkerSize', 10);
hold on
plot(saveTime, saveC(:, 2), '-*', 'MarkerSize', 10);
hold on
plot(saveTime, saveC(:, 3), '-*', 'MarkerSize', 10);
set(gca, 'fontsize', 15);
legend('A', 'B', 'C');
hold on
plot(saveTime, repmat(1.4056, numel(saveTime), 1))
hold on
plot(saveTime, repmat(0.7114, numel(saveTime), 1))
%title('Concentrations');
xlabel('t');

%Make spatial plots
evaluateInterface(cellGrid, Xi, distFunctions, false);
toc

levelSet = dV;
levelSet(Xi == 1) = -levelSet(Xi == 1);
levelSet(Xi == 2) = max(levelSet(Xi == 2), 10^(-4));
levelSet(Xi == 3) = max(levelSet(Xi == 3), 10^(-4));
gridHyPHMBasis = Grid(coord, cellGrid.triangles);
[gridHyPHM, levelSetNew] = localMeshRefinementTripel(gridHyPHMBasis, ...
    levelSet);

threshold = 1 * 10^(-15);
temp = c(1) * ones(gridHyPHM.numT, 1);
markBasis = sum(reshape(Xi(gridHyPHMBasis.V0T(:, :), 1), gridHyPHMBasis.numT, 3) == 3, 2) > 0.5;
markSolid = all(levelSetNew(gridHyPHM.V0T) > 0, 2);
Interp = scatteredInterpolant(gridHyPHMBasis.baryT(:, 1), gridHyPHMBasis.baryT(:, 2), double(markBasis), 'nearest');
markRefine = logical(Interp(gridHyPHM.baryT(:, 1), gridHyPHM.baryT(:, 2)));
temp(markRefine & markSolid) = threshold;

threshold = 0;
markBasis = sum(reshape(Xi(gridHyPHMBasis.V0T(:, :), 1), gridHyPHMBasis.numT, 3) == 2, 2) > 0.5;
markSolid = all(levelSetNew(gridHyPHM.V0T) > 0, 2);
Interp = scatteredInterpolant(gridHyPHMBasis.baryT(:, 1), gridHyPHMBasis.baryT(:, 2), double(markBasis), 'nearest');
markRefine = logical(Interp(gridHyPHM.baryT(:, 1), gridHyPHM.baryT(:, 2)));
temp(markRefine & markSolid) = threshold;

Aconcentration = Variable(gridHyPHM, Stepper([0, 1]), ...
    'A', 'P0');
Aconcentration.setdata(temp);
%Aconcentration.visualize;

%% Symbolic computation (equilibrium point, linear stability analysis)

syms a b c d e
eqns = [a / b - 1 == 0; c * b - 1 == 0; a * (1 - d - e) - 20 * e - 2 * (1 - 0.04 * pi) + 0.4 * pi == 0; b * (1 - d - e) + 20 * e + d * 4 - (1 - 0.04 * pi + 24 * 0.02 * pi) == 0; c * (1 - d - e) + 4 * d - (1 - 0.04 * pi + 4 * 0.02 * pi) == 0];
S = solve(eqns, [a, b, c, d, e]);
answer = [S.a(2, 1), S.b(2, 1), S.c(2, 1), S.d(2, 1), S.e(2, 1)];
s = answer;

f = @(x) [x(1) / x(2) - 1; x(3) * x(2) - 1; x(1) * (1 - x(4) - x(5)) - 20 * x(5) - 2 * (1 - 0.04 * pi) + 0.4 * pi; x(2) * (1 - x(4) - x(5)) + 20 * x(5) + x(4) * 4 - (1 - 0.04 * pi + 24 * 0.02 * pi); x(3) * (1 - x(4) - x(5)) + 4 * x(4) - (1 - 0.04 * pi + 4 * 0.02 * pi)];
solution = fsolve(f, ones(5, 1), options);
[saveC(end, :)'; saveSolidArea(end) - saveDArea(end); saveDArea(end)]

f = @(x) [s(5) * (x(2) / x(1) - 1) * (20 + x(1)) + s(4) * (x(2) * x(3) - 1) * x(1); ...
    -s(4) * (x(2) * x(3) - 1) * (4 - x(2)) - s(5) * (x(2) / x(1) - 1) * (20 - x(2)); ...
    -s(4) * (x(2) * x(3) - 1) * (4 - x(3)) + s(5) * (x(2) / x(1) - 1) * x(3)];

A = -[-s(5) * s(2) / s(1)^2 * (-20 - s(1)), s(5) / s(1) * (-20 - s(1)) - s(4) * s(3) * s(1), -s(4) * s(2) * s(1); ...
    -s(5) * s(2) / s(1)^2 * (20 - s(2)), s(5) / s(1) * (20 - s(2)) + s(4) * s(3) * (4 - s(2)), s(4) * s(2) * (4 - s(2)); ...
    s(5) * s(2) / s(1)^2 * s(3), -s(5) / s(1) * s(3) + s(4) * s(3) * (4 - s(3)), s(4) * s(2) * (4 - s(3))];
