% Simple test script for 2 mineral dissolution without VIIM at the
% porescale

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% General setting variables
clear initialLevelSetDataCells
global Solver
Solver = 'StandardDirect';
global EPS;
EPS = eps;


lengthXAxis = 0.1; % [ dm ]
lengthYAxis = 0.1; % [ dm ]
numPartitionsMicroscale = 64; % Number of partitions
macroscaleStepSize = lengthXAxis / 5; % [ dm ]


tic; % Preprocessing

% Physical parameters
dimension = 2;
spaceScaleFactor = 1; % spaceScaleFactor [length of Y] = 1 [cm]

% Parameters in dm and s
diffusionCoefficient = 0.2 * 1e-4 * 1e-2; % [ dm^2 s^(-1) ]

molarVolumeC = 39.63 * 1e-3;
molarVolumeD = 64.12 * 1e-3; % [ dm^3 mol^(-1) ]

rateCoefficientHydrogen = 0.89 * 1e-2; % [ mol dm^(-2) s^(-1) ]

rateCoefficientTSTC = 6.6e-7 * 1e-2; % [ mol dm^(-2) s^(-1) ]
rateCoefficientTSTD = 4.5e-4 * 1e-2; % [ mol dm^(-2) s^(-1) ]

equilibriumRateConstantC = 10^(-8.234);
equilibriumRateConstantD = 10^(-16.5);
%inletVelocity =0.1* 1e-1; % [ dm s^(-1) ]; Pe = 5000

%inletVelocity = 0.001 * 1e-1; % [ dm s^(-1) ]; Pe = 50
inletVelocity = 0.01 * 1e-1; % [ dm s^(-1) ]; Pe = 500
initialHydrogenConcentration = 10^(-5); % [ mol dm^(-3) ]

dissolutionReactionRateC = @(cHydrogen, cCalcium, cCarbonate) ...
    (rateCoefficientHydrogen * cHydrogen + rateCoefficientTSTC) ...
    * (1 - cCalcium * cCarbonate / equilibriumRateConstantC);

dissolutionReactionRateD = @(cHydrogen, cCalcium, cCarbonate, cMagnesium) ...
    (cHydrogen^0.5 * rateCoefficientTSTD) ...
    * (1 - cCalcium * cCarbonate^2 * cMagnesium / equilibriumRateConstantD);


numberOfSlices = 1;
disp(['numSlices = ', num2str(numberOfSlices)]);
pecletNumber = (inletVelocity * lengthXAxis) / diffusionCoefficient;
fprintf(['Peclet number: ', num2str(pecletNumber), '\n']);

% dt_macro = 1 * 3600; % Pe = 5000
% initialMacroscaleTimeStepSize = 160; % Pe = 50
initialMacroscaleTimeStepSize = 10;
% dt_micro = dt_macro / 100; % 0;

endTime = 1 * 3600 * 600; % [h] Pe = 5000
% endTime = 10 * 60 * 60; % Pe = 50
% endTime = 200 * dt_macro;

%% Computation of time steps in simulation
%timeStepperType = 'exp';
timeStepperType = 'expmax';

switch (timeStepperType)
    case 'linear'
        if (mod(endTime, initialMacroscaleTimeStepSize) < EPS)
            timeSteps = 0:initialMacroscaleTimeStepSize:endTime;
        else
            timeSteps = [0:initialMacroscaleTimeStepSize:endTime, endTime];
        end
        timeStepSizeFactor = 1;
    case 'exp'
        %timeStepSizeFactor = sqrt(2);
        timeStepSizeFactor = 2;
        numTimeSteps = floor(log(endTime / initialMacroscaleTimeStepSize) ...
            /log(timeStepSizeFactor));
        timeSteps = [0; timeStepSizeFactor.^(0:numTimeSteps)'] ...
            * initialMacroscaleTimeStepSize;
        if (timeSteps(end) < endTime)
            timeSteps(end+1) = endTime;
        end

    case 'expmax'
        maximalStep = 0.8 * 10^5;
        timeStepSizeFactor = 2;
        timeSteps = [0, initialMacroscaleTimeStepSize];
        while timeSteps(end) < endTime
            timeSteps = [timeSteps, min(timeSteps(end) * timeStepSizeFactor, timeSteps(end) + maximalStep)];
        end
        timeSteps = [timeSteps(1:(end -1)), endTime];
        size(timeSteps)

    otherwise
        error('Time stepper type not implemented.');
end
numTimeSlices = numel(timeSteps);
levelSetEvolutionTime = NaN(numTimeSlices, 1);
cellProblemTime = NaN(numTimeSlices, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization of level set method variables

microscaleGrid = FoldedCartesianGrid(dimension, ...
    [0, lengthXAxis, 0, lengthYAxis], ...
    numPartitionsMicroscale*[round(lengthXAxis / lengthYAxis), 1]);
coord = microscaleGrid.coordinates;
coordCell = mat2cell(coord, ones(1, microscaleGrid.nodes), dimension);


initialLevelSetFunc = @(x) 0.03 - norm(x-[0.05, 0.05]);
coordCell = mat2cell(coord, ones(1, microscaleGrid.nodes), dimension);
initialLevelSetDataCells{1}(:, 1) = cellfun(initialLevelSetFunc, coordCell);

[a, b] = meshgrid(0:lengthYAxis/numPartitionsMicroscale:lengthXAxis, 0:lengthYAxis/numPartitionsMicroscale:lengthYAxis);
contour(a, b, reshape(initialLevelSetDataCells{1}, numPartitionsMicroscale + 1, numPartitionsMicroscale + 1)', [0, 1])
axis equal

interfaceNormalVelocityC = @(cHydrogen, cCalcium, cCarbonate) ...
    spaceScaleFactor * molarVolumeC * dissolutionReactionRateC(cHydrogen, ...
    cCalcium, cCarbonate);

interfaceNormalVelocityD = @(cHydrogen, cCalcium, cCarbonate, cMagnesium) ...
    spaceScaleFactor * molarVolumeD * dissolutionReactionRateD(cHydrogen, ...
    cCalcium, cCarbonate, cMagnesium);

currentTime = 0;

currentLevelSetDataCells = cell(numberOfSlices, 1);
[currentLevelSetDataCells{:}] = deal(initialLevelSetDataCells{1});
oldLevelSetDataCells = currentLevelSetDataCells;
% TODO Viel zu viele Spalten in levelSetData
% prev: ( ceil( endTime / dt_macro ) + 1 ) columns
levelSetData = NaN(numel(initialLevelSetDataCells{1}), numTimeSlices);
levelSetData(:, 1) = currentLevelSetDataCells{1};
levelSet = cell(numberOfSlices, 1);
[levelSet{:}] = deal(levelSetData);
clear levelSetData;

gridHyPHM = Grid(coord, microscaleGrid.triangles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i = 1:gridHyPHM.numE
    if (all(gridHyPHM.coordV(gridHyPHM.V0E(i, :), 1) < eps))
        gridHyPHM.idE(i) = 4;
    end
    if (all(gridHyPHM.coordV(gridHyPHM.V0E(i, :), 1) > lengthXAxis - eps))
        gridHyPHM.idE(i) = 2;
    end
    if (all(gridHyPHM.coordV(gridHyPHM.V0E(i, :), 2) < eps))
        gridHyPHM.idE(i) = 1;
    end
    if (all(gridHyPHM.coordV(gridHyPHM.V0E(i, :), 2) > lengthYAxis - eps))
        gridHyPHM.idE(i) = 3;
    end
end

macroCoordCell = mat2cell(gridHyPHM.baryT, ones(gridHyPHM.numT, 1), 2);

% timeSteps(end) = []; % At the last time step, no computations are neccessary.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculation of initial effective parameters

diffusionTensors = cell(numberOfSlices, 1);
[diffusionTensors{:}] = deal(NaN(4, numTimeSlices));
porosities = cell(numberOfSlices, 1);
[porosities{1:end}] = deal(NaN(numTimeSlices, 1));
clear numTimeSlices;

surfaceArea = cell(numberOfSlices, 1);
PorositiyEstimate = sum(levelSet{1}(:, 1) < -eps) / numel(levelSet{1}(:, 1))

for i = 1:numberOfSlices

    [cellProblemSystemMatrix, rhs, isDoF, triangleVolumes, triangleSurfaces] ...
        = assembleCellProblem(microscaleGrid, levelSet{i}(:, 1));

    surfaceArea{i}(1) = sum(triangleSurfaces);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization of HyPHM variables (cont.)

flowStepper = Stepper(0:1);
transportStepper = Stepper(timeSteps);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Stokes
%Calculate Pressure distribution and Stokes velocity
disp(' ');
disp('Initializing Stokes problem...');

phead = Variable(gridHyPHM, flowStepper, 'pressure head', 'P1');
StokesVelocity = Variable(gridHyPHM, flowStepper, 'Stokes Velocity', 'P2P2');
phead.setdata(0, @(t, x) 0.0);
StokesVelocity.setdata(0, @(t, x) 0.0);

StokesL = StokesLEVEL(gridHyPHM, flowStepper, 'Stokes problem');
StokesL.L.setdata(levelSet{1}(:, 1))
StokesL.id2D = {4, 3, 1};
StokesL.uD.setdata(@(t, x) inletVelocity*(x(1) < EPS)*[1; 0]);
StokesL.F.setdata(@(t, x) 0);
StokesL.U = StokesVelocity;
StokesL.P = phead;

flowStepper.next;
StokesL.computeLevel('s');
flowStepper.prev;
%   StokesL.U.visualize();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


flow = Variable(gridHyPHM, flowStepper, 'Flow', 'RT0');
helper = StokesL.U.getdata(1);
flow.setdata(helper(gridHyPHM.numV + 1:end, 1).*gridHyPHM.nuE(:, 1)+helper(gridHyPHM.numV + 1:end, 2).*gridHyPHM.nuE(:, 2));


% Boundary IDs: 1 = down, 2 = right, 3 = up, 4 = left

disp(' ');
disp('Initializing calcium transport...');


calciumConcentration = Variable(gridHyPHM, transportStepper, ...
    'Ca^(2+)', 'P0');
calciumConcentration.setdata(0, 0*ones(gridHyPHM.numT, 1));
% calciumConcentration.setdata( 0, @(t,x) sqrt( equilibriumRateConstant ) );
calciumTransport = TransportLEVEL(gridHyPHM, transportStepper, 'Ca^(2+) Transport');
% calciumTransport.id2D = {2};
calciumTransport.D.setdata(diffusionCoefficient*eye(2));
calciumTransport.id2N = {1, 2, 3};
calciumTransport.id2F = {4};
calciumTransport.U = calciumConcentration;
% calciumTransport.uD = Variable( g_HyPHM, st, 'uD', 'P0E' );
% calciumTransport.uD.setdata( zeros( g_HyPHM.numE, 1 ) );
calciumTransport.gF = Variable(gridHyPHM, transportStepper, 'gF', 'P0E');
calciumTransport.gF.setdata(zeros(gridHyPHM.numE, 1));
calciumTransport.A = Variable(gridHyPHM, transportStepper, 'Porosity', 'P0');
calciumTransport.A.setdata(0, ones(gridHyPHM.numT, 1));
% calciumTransport.C.setdata( P2P2.P2P2toRT0slice( g_HyPHM, flow.getdata(1) ) );
calciumTransport.C.setdata(0, flow.getdata(1));
calciumTransport.isUpwind = 'exp';

disp(' ');
disp('Initializing carbonate transport...');

carbonateConcentration = Variable(gridHyPHM, transportStepper, ...
    'CO_3^(2-)', 'P0');
carbonateConcentration.setdata(0, 0*ones(gridHyPHM.numT, 1));
% carbonateConcentration.setdata( 0, @(t,x) sqrt( equilibriumRateConstant ) );
carbonateTransport = TransportLEVEL(gridHyPHM, transportStepper, 'CO_3^(2-) Transport');
% carbonateTransport.id2D = {2};
carbonateTransport.D.setdata(diffusionCoefficient*eye(2));
carbonateTransport.id2N = {1, 2, 3};
carbonateTransport.id2F = {4};
carbonateTransport.U = carbonateConcentration;
% carbonateTransport.uD = Variable( g_HyPHM, st, 'uD', 'P0E' );
% carbonateTransport.uD.setdata( zeros( g_HyPHM.numE, 1 ) );
carbonateTransport.gF = Variable(gridHyPHM, transportStepper, 'gF', 'P0E');
carbonateTransport.gF.setdata(zeros(gridHyPHM.numE, 1));
carbonateTransport.A = Variable(gridHyPHM, transportStepper, 'Porosity', 'P0');
carbonateTransport.A.setdata(0, ones(gridHyPHM.numT, 1));

% carbonateTransport.C.setdata( P2P2.P2P2toRT0slice( g_HyPHM, flow.getdata(1) ) );
carbonateTransport.C.setdata(flow.getdata(1));
carbonateTransport.isUpwind = 'exp';

disp(' ');
disp('Initializing magnesium transport...');

magnesiumConcentration = Variable(gridHyPHM, transportStepper, ...
    'Mg^(2+)', 'P0');
magnesiumConcentration.setdata(0, 0*ones(gridHyPHM.numT, 1));
% magnesiumConcentration.setdata( 0, @(t,x) sqrt( equilibriumRateConstant ) );
magnesiumTransport = TransportLEVEL(gridHyPHM, transportStepper, 'MG^(2+) Transport');
% magnesiumTransport.id2D = {2};
magnesiumTransport.D.setdata(diffusionCoefficient*eye(2));
magnesiumTransport.id2N = {1, 2, 3};
magnesiumTransport.id2F = {4};
magnesiumTransport.U = magnesiumConcentration;
% magnesiumTransport.uD = Variable( g_HyPHM, st, 'uD', 'P0E' );
% magnesiumTransport.uD.setdata( zeros( g_HyPHM.numE, 1 ) );
magnesiumTransport.gF = Variable(gridHyPHM, transportStepper, 'gF', 'P0E');
magnesiumTransport.gF.setdata(zeros(gridHyPHM.numE, 1));
magnesiumTransport.A = Variable(gridHyPHM, transportStepper, 'Porosity', 'P0');
magnesiumTransport.A.setdata(0, ones(gridHyPHM.numT, 1));
% magnesiumTransport.C.setdata( P2P2.P2P2toRT0slice( g_HyPHM, flow.getdata(1) ) );
magnesiumTransport.C.setdata(0, flow.getdata(1));
magnesiumTransport.isUpwind = 'exp';
disp(' ');
disp('Initializing hydrogen transport...');

hydrogenConcentration = Variable(gridHyPHM, transportStepper, ...
    'H^+_pH', 'P0');
hydrogenConcentration.setdata(0, initialHydrogenConcentration*ones(gridHyPHM.numT, 1));
hydrogenTransport = TransportLEVEL(gridHyPHM, transportStepper, 'H^+ Transport');
hydrogenTransport.id2N = {1, 2, 3};
hydrogenTransport.id2F = {4};
% hydrogenTransport.id2D = {4};
hydrogenTransport.U = hydrogenConcentration;
hydrogenTransport.D.setdata(diffusionCoefficient*eye(2));
hydrogenTransport.gF.setdata( ...
    @(t, x) -initialHydrogenConcentration*inletVelocity*(x(1) < EPS));

% hydrogenTransport.gF.setdata( @(t,x) -1e-7 );
% hydrogenTransport.uD.setdata( @(t,x) initialHydrogenConcentration );
hydrogenTransport.A.setdata(0, ones(gridHyPHM.numT, 1));
hydrogenTransport.C.setdata(0, flow.getdata(1));
hydrogenTransport.isUpwind = 'exp';

speciesCells = {calciumTransport; carbonateTransport; hydrogenTransport; magnesiumTransport};

calciumDataOld = calciumConcentration.getdata(0);
carbonateDataOld = carbonateConcentration.getdata(0);
magneisumDataOld = magnesiumConcentration.getdata(0);
hydrogenDataOld = hydrogenConcentration.getdata(0);

preprocessingTime = toc; % Preprocessing

disp(' ');
disp(['Preprocessing done in ', num2str(preprocessingTime), ' seconds.']);


hydrogenTransport.L.setdata(0, [levelSet{1}(:, 1)]);
magnesiumTransport.L.setdata(0, [levelSet{1}(:, 1)]);
carbonateTransport.L.setdata(0, [levelSet{1}(:, 1)]);
calciumTransport.L.setdata(0, [levelSet{1}(:, 1)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mapping Vertices --> traingles

VertexTriMatrix = cell(1, gridHyPHM.numV);
for i = 1:gridHyPHM.numT
    for j = 1:3
        VertexTriMatrix{1, gridHyPHM.V0T(i, j)} = [VertexTriMatrix{1, gridHyPHM.V0T(i, j)}, i];
    end
end

%% Time iteration

%while ( t_ <= endTime + EPS )
while transportStepper.next


    timeIterationStep = transportStepper.curstep;
    macroscaleTimeStepSize = transportStepper.curtau;
    disp('----------------------------------------');
    disp(['Time step ', num2str(timeIterationStep)]);
    disp(['  Current time = ', ...
        num2str(currentTime)]);
    disp(['  Current time step size = ', ...
        num2str(macroscaleTimeStepSize)]);

    if (currentTime + macroscaleTimeStepSize > endTime - EPS)
        macroscaleTimeStepSize = endTime - currentTime;
        if (abs(macroscaleTimeStepSize) < EPS)
            break;
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   inner corrector
    speciesCells = {calciumTransport; carbonateTransport; hydrogenTransport; magnesiumTransport};
    [out, Corrector] = Continuation(speciesCells, timeIterationStep, levelSet{1}(:, timeIterationStep), 20);
    calciumTransport = out{1};
    carbonateTransport = out{2};
    hydrogenTransport = out{3};
    magnesiumTransport = out{4};


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Level set evolution step

    disp('Evolution of level set ...');
    tic;
    coordCell = mat2cell(microscaleGrid.coordinates, ...
        ones(1, microscaleGrid.nodes), dimension);
    timeCell = num2cell(currentTime*ones(microscaleGrid.nodes, 1));


    calciumData = calciumTransport.U.getdata(timeIterationStep-1);
    carbonateData = carbonateTransport.U.getdata(timeIterationStep-1);
    magnesiumData = magnesiumTransport.U.getdata(timeIterationStep-1);
    hydrogenData = max(hydrogenTransport.U.getdata(timeIterationStep - 1), eps);


    %Convert Concentrations from T --> V
    calciumDataV = zeros(gridHyPHM.numV, 1);
    magnesiumDataV = zeros(gridHyPHM.numV, 1);
    hydrogenDataV = zeros(gridHyPHM.numV, 1);
    carbonateDataV = zeros(gridHyPHM.numV, 1);

    for i = 1:gridHyPHM.numV
        %[row,~] = find(gridHyPHM.V0T(:,:)==i);
        row = VertexTriMatrix{1, i};
        calciumDataV(i) = max(calciumData(row));
        magnesiumDataV(i) = max(magnesiumData(row));
        hydrogenDataV(i) = max(hydrogenData(row));
        carbonateDataV(i) = max(carbonateData(row));
    end

    normalSpeedC = arrayfun(interfaceNormalVelocityC, ...
        hydrogenDataV, ...
        calciumDataV, ...
        carbonateDataV);

    normalSpeedD = arrayfun(interfaceNormalVelocityD, ...
        hydrogenDataV, ...
        calciumDataV, ...
        carbonateDataV, ...
        magnesiumDataV);


    normalSpeed = normalSpeedC .* (gridHyPHM.coordV(:, 2) > lengthYAxis / 2) + normalSpeedD .* (gridHyPHM.coordV(:, 2) <= lengthYAxis / 2);

    InterfacePoints = getInterfacePoints(levelSet{1}(:, timeIterationStep), numPartitionsMicroscale*lengthXAxis);
    normalSpeednew = zeros(numel(levelSet{1}(:, timeIterationStep)), numel(InterfacePoints));
    for i = 1:numel(InterfacePoints)
        center = InterfacePoints(i);
        neighbors = getClosePoints(center, 0.01, gridHyPHM);
        normalSpeednew(neighbors, i) = normalSpeed(center);
    end
    normalSpeed = sum(normalSpeednew, 2) ./ max(sum(normalSpeednew > 0, 2), 1);
    imagesc(reshape(normalSpeed, 65, 65)')
    %surface(a,b,reshape(normalSpeed,71,71)')

    for i = 1:numberOfSlices


        CFL = 1 / 4 * 1 / 10 * 1 / numPartitionsMicroscale / max(normalSpeed);
        %microscaleTimeStepSize = transportStepper.curtau / numMicroscaleTimeSteps;
        oldMicroscaleTime = currentTime;
        for j = 1:ceil(transportStepper.curtau/CFL)
            microscaleTimeStepSize = min(currentTime+transportStepper.curtau-oldMicroscaleTime, CFL);
            %             disp( ['    Level set substep ', num2str(j) ] );

            newMicroscaleTime = oldMicroscaleTime + microscaleTimeStepSize;
            % [] argument is unused in method (needed for implicit methods)
            currentLevelSetDataCells{i} = levelSetEquationTimeStep( ...
                newMicroscaleTime, ...
                oldMicroscaleTime, oldLevelSetDataCells{i}, microscaleGrid, normalSpeed);
            oldMicroscaleTime = newMicroscaleTime;
            oldLevelSetDataCells{i} = currentLevelSetDataCells{i};

        end
        j
        levelSet{i}(:, timeIterationStep + 1) = currentLevelSetDataCells{i};
    end


    currentTime = currentTime + macroscaleTimeStepSize;

    levelSetEvolutionTime(timeIterationStep) = toc;
    disp(['    ... done in ', ...
        num2str(levelSetEvolutionTime(timeIterationStep)), ' seconds.']);

    hydrogenTransport.L.setdata(timeIterationStep, levelSet{1}(:, timeIterationStep + 1));
    magnesiumTransport.L.setdata(timeIterationStep, levelSet{1}(:, timeIterationStep + 1));
    carbonateTransport.L.setdata(timeIterationStep, levelSet{1}(:, timeIterationStep + 1));
    calciumTransport.L.setdata(timeIterationStep, levelSet{1}(:, timeIterationStep + 1));


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Stokes
    %Calculate Pressure distribution and Stokes velocity
    disp(' ');
    disp('Initializing Stokes problem...');


    StokesL.L.setdata(levelSet{1}(:, timeIterationStep + 1))

    flowStepper.next;
    StokesL.computeLevel('s');

    helper = StokesL.U.getdata(1);
    flow.setdata(helper(gridHyPHM.numV + 1:end, 1).*gridHyPHM.nuE(:, 1)+helper(gridHyPHM.numV + 1:end, 2).*gridHyPHM.nuE(:, 2));
    flowStepper.prev;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Macroscopic transport step


    Levels = hydrogenTransport.L.getdata(timeIterationStep);

    SorceScale = zeros(gridHyPHM.numT, 1);

    for i = 1:3
        temp1 = find(Levels(gridHyPHM.V0E(gridHyPHM.E0T(:, i), 1)) > -eps & Levels(gridHyPHM.V0E(gridHyPHM.E0T(:, i), 2)) > -eps);
        SorceScale(temp1) = SorceScale(temp1) + gridHyPHM.areaE(gridHyPHM.E0T(temp1, i)) ./ gridHyPHM.areaT(temp1);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nonlinearFunc = cell(4, 1);
    nonlinearFunc{1} = @(x, y, z, w) (1 ...
        .* (sqrt(max(z, eps)) * rateCoefficientTSTD) ...
        .* (1 - x .* y .* y .* w ./ equilibriumRateConstantD) .* (gridHyPHM.baryT(:, 2) <= lengthYAxis / 2) ...
        +(rateCoefficientHydrogen * z + rateCoefficientTSTC) .* ...
        (1 - x .* y ./ equilibriumRateConstantC) .* (gridHyPHM.baryT(:, 2) > lengthYAxis / 2)) .* SorceScale;

    nonlinearFunc{2} = @(x, y, z, w) (2 ...
        .* (sqrt(max(z, eps)) * rateCoefficientTSTD) ...
        .* (1 - x .* y .* y .* w ./ equilibriumRateConstantD) .* (gridHyPHM.baryT(:, 2) <= lengthYAxis / 2) ...
        +(rateCoefficientHydrogen * z + rateCoefficientTSTC) .* ...
        (1 - x .* y ./ equilibriumRateConstantC) .* (gridHyPHM.baryT(:, 2) > lengthYAxis / 2)) .* SorceScale;
    nonlinearFunc{3} = @(x, y, z, w) -1 ...
        .* rateCoefficientHydrogen .* z ...
        .* (1 - x .* y ./ equilibriumRateConstantC) .* (gridHyPHM.baryT(:, 2) > lengthYAxis / 2) .* SorceScale;
    nonlinearFunc{4} = @(x, y, z, w) 1 ...
        .* (sqrt(max(z, eps)) * rateCoefficientTSTD) ...
        .* (1 - x .* y .* y .* w ./ equilibriumRateConstantD) .* SorceScale .* (gridHyPHM.baryT(:, 2) <= lengthYAxis / 2);

    nonlinearJacFunc = cell(4, 4);
    nonlinearJacFunc{1, 1} = @(x, y, z, w) (1 ...
        .* (sqrt(max(z, eps)) * rateCoefficientTSTD) ...
        .* (-y .* y .* w ./ equilibriumRateConstantD) .* (gridHyPHM.baryT(:, 2) <= lengthYAxis / 2) ...
        +(rateCoefficientHydrogen * z + rateCoefficientTSTC) .* ...
        (-y ./ equilibriumRateConstantC) .* (gridHyPHM.baryT(:, 2) > lengthYAxis / 2)) .* SorceScale;
    nonlinearJacFunc{1, 2} = @(x, y, z, w) (1 ...
        .* (sqrt(max(z, eps)) * rateCoefficientTSTD) ...
        .* (-2 * x .* y .* w ./ equilibriumRateConstantD) .* (gridHyPHM.baryT(:, 2) <= lengthYAxis / 2) ...
        +(rateCoefficientHydrogen * z + rateCoefficientTSTC) .* ...
        (-x ./ equilibriumRateConstantC) .* (gridHyPHM.baryT(:, 2) > lengthYAxis / 2)) .* SorceScale;
    nonlinearJacFunc{1, 3} = @(x, y, z, w) (1 ...
        .* (0.5 ./ sqrt(max(z, eps)) * rateCoefficientTSTD) ...
        .* (1 - x .* y .* y .* w ./ equilibriumRateConstantD) .* (gridHyPHM.baryT(:, 2) <= lengthYAxis / 2) ...
        +(rateCoefficientHydrogen) .* ...
        (1 - x .* y ./ equilibriumRateConstantC) .* (gridHyPHM.baryT(:, 2) > lengthYAxis / 2)) .* SorceScale;
    nonlinearJacFunc{1, 4} = @(x, y, z, w) 1 ...
        .* (sqrt(max(z, eps)) * rateCoefficientTSTD) ...
        .* (-x .* y .* y ./ equilibriumRateConstantD) .* (gridHyPHM.baryT(:, 2) <= lengthYAxis / 2) .* SorceScale;

    nonlinearJacFunc{2, 1} = @(x, y, z, w) (2 ...
        .* (sqrt(max(z, eps)) * rateCoefficientTSTD) ...
        .* (-y .* y .* w ./ equilibriumRateConstantD) .* (gridHyPHM.baryT(:, 2) <= lengthYAxis / 2) ...
        +(rateCoefficientHydrogen * z + rateCoefficientTSTC) .* ...
        (-y ./ equilibriumRateConstantC) .* (gridHyPHM.baryT(:, 2) > lengthYAxis / 2)) .* SorceScale;
    nonlinearJacFunc{2, 2} = @(x, y, z, w) (2 ...
        .* (sqrt(max(z, eps)) * rateCoefficientTSTD) ...
        .* (-2 * x .* y .* w ./ equilibriumRateConstantD) .* (gridHyPHM.baryT(:, 2) <= lengthYAxis / 2) ...
        +(rateCoefficientHydrogen * z + rateCoefficientTSTC) .* ...
        (-x ./ equilibriumRateConstantC) .* (gridHyPHM.baryT(:, 2) > lengthYAxis / 2)) .* SorceScale;
    nonlinearJacFunc{2, 3} = @(x, y, z, w) (2 ...
        .* (0.5 ./ sqrt(max(z, eps)) * rateCoefficientTSTD) ...
        .* (1 - x .* y .* y .* w ./ equilibriumRateConstantD) .* (gridHyPHM.baryT(:, 2) <= lengthYAxis / 2) ...
        +(rateCoefficientHydrogen) .* ...
        (1 - x .* y ./ equilibriumRateConstantC) .* (gridHyPHM.baryT(:, 2) > lengthYAxis / 2)) .* SorceScale;
    nonlinearJacFunc{2, 4} = @(x, y, z, w) 2 ...
        .* (sqrt(max(z, eps)) * rateCoefficientTSTD) ...
        .* (-x .* y .* y ./ equilibriumRateConstantD) .* (gridHyPHM.baryT(:, 2) <= lengthYAxis / 2) .* SorceScale;

    nonlinearJacFunc{3, 1} = @(x, y, z, w) -1 ...
        .* rateCoefficientHydrogen .* z ...
        .* (-y ./ equilibriumRateConstantC) .* (gridHyPHM.baryT(:, 2) > lengthYAxis / 2) .* SorceScale;
    nonlinearJacFunc{3, 2} = @(x, y, z, w) -1 ...
        .* rateCoefficientHydrogen .* z ...
        .* (-x ./ equilibriumRateConstantC) .* (gridHyPHM.baryT(:, 2) > lengthYAxis / 2) .* SorceScale;
    nonlinearJacFunc{3, 3} = @(x, y, z, w) -1 ...
        .* rateCoefficientHydrogen ...
        .* (1 - x .* y ./ equilibriumRateConstantC) .* (gridHyPHM.baryT(:, 2) > lengthYAxis / 2) .* SorceScale;
    nonlinearJacFunc{3, 4} = @(x, y, z, w) 0;

    nonlinearJacFunc{4, 1} = @(x, y, z, w) 1 ...
        .* (sqrt(max(z, eps)) * rateCoefficientTSTD) ...
        .* (-y .* y .* w ./ equilibriumRateConstantD) .* SorceScale .* (gridHyPHM.baryT(:, 2) <= lengthYAxis / 2);
    nonlinearJacFunc{4, 2} = @(x, y, z, w) 1 ...
        .* (sqrt(max(z, eps)) * rateCoefficientTSTD) ...
        .* (-2 * x .* y .* w ./ equilibriumRateConstantD) .* SorceScale .* (gridHyPHM.baryT(:, 2) <= lengthYAxis / 2);
    nonlinearJacFunc{4, 3} = @(x, y, z, w) 1 ...
        .* (0.5 ./ sqrt(max(z, eps)) * rateCoefficientTSTD) ...
        .* (1 - x .* y .* y .* w ./ equilibriumRateConstantD) .* SorceScale .* (gridHyPHM.baryT(:, 2) <= lengthYAxis / 2);
    nonlinearJacFunc{4, 4} = @(x, y, z, w) 1 ...
        .* (sqrt(max(z, eps)) * rateCoefficientTSTD) ...
        .* (-x .* y .* y ./ equilibriumRateConstantD) .* SorceScale .* (gridHyPHM.baryT(:, 2) <= lengthYAxis / 2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    speciesCells = {calciumTransport; carbonateTransport; hydrogenTransport; magnesiumTransport};

    calciumTransport.A.setdata(timeIterationStep, ones(gridHyPHM.numT, 1));
    calciumTransport.C.setdata(timeIterationStep, flow.getdata(1));
    carbonateTransport.A.setdata(timeIterationStep, ones(gridHyPHM.numT, 1));
    carbonateTransport.C.setdata(timeIterationStep, flow.getdata(1));
    hydrogenTransport.A.setdata(timeIterationStep, ones(gridHyPHM.numT, 1));
    hydrogenTransport.C.setdata(timeIterationStep, flow.getdata(1));
    magnesiumTransport.A.setdata(timeIterationStep, ones(gridHyPHM.numT, 1));
    magnesiumTransport.C.setdata(timeIterationStep, flow.getdata(1));
    tic
    newtonIteration(speciesCells, nonlinearFunc, nonlinearJacFunc, 10);
    toc


    calciumTransport.U.setdata(timeIterationStep-1, calciumTransport.U.getdata(timeIterationStep - 1)-Corrector(:, 1));
    carbonateTransport.U.setdata(timeIterationStep-1, carbonateTransport.U.getdata(timeIterationStep - 1)-Corrector(:, 2));
    hydrogenTransport.U.setdata(timeIterationStep-1, max(hydrogenTransport.U.getdata(timeIterationStep - 1) - Corrector(:, 3), eps));
    magnesiumTransport.U.setdata(timeIterationStep-1, magnesiumTransport.U.getdata(timeIterationStep - 1)-Corrector(:, 4));


end % while

calciumConcentration.visualize;
% carbonateConcentration.visualize;
%hydrogenTransport.U.visualize;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Helper functions


function diff = diffusionHelperFun(t, x, diffusionCells, numberOfSlices, timeStep)

diff = zeros(2);
for i = 1:numberOfSlices
    diff = diff + reshape(diffusionCells{1}(:, timeStep), 2, 2);
end

%     diff = reshape( diffusionCells{ findSlice( x, numberOfSlices ) }( :, ...
%         timeStep ), 2, 2 );

end


function out = getInterfacePoints(Levels, dist)

out = find((Levels < 0).*(Levels > -1 / (dist)));

end

function out = getClosePoints(Center, dist, grid)
out = find(sqrt((grid.coordV(:, 1)-grid.coordV(Center, 1)).^2 + (grid.coordV(:, 2) - grid.coordV(Center, 2)).^2) < dist);
end
