% Perform micro-macro simulation with dolomite grains.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% General setting variables
global Solver
Solver = 'StandardIterative';
all = tic;
global EPS;
EPS = eps;
global timeNewton;
global timeLevelSet;
global timeDiff;
global timePerm;
global timeDracy;
global numPartitionsMicroscale;
global macroscaleStepSize;
global PermeabilityTensorTest;
global CalciumConcentrationTest;
global GridTest;

timeNewton = 0;
timeLevelSet = 0;
timeDiff = 0;
timePerm = 0;
timeDracy = 0;
numProc = 1;
lengthXAxis = 1 * 0.1; % [ dm ]
lengthYAxis = 0.5 * 0.1; % [ dm ]
numPartitionsMicroscale = 64; % Number of partitions in each direction
macroscaleStepSize = lengthXAxis / (27 * 4); % [ dm ]
meanGrainSize = 10^(-4); % [ dm ]

tic; % Preprocessing

% Physical parameters
dimension = 2;
spaceScaleFactor = 10 * sqrt(13*27); % spaceScaleFactor [length of Y] = 1 [cm]

% Parameters in dm and s
dynamicViscosity = 1100e-7; % [ kg dm^(-1) s^(-1) ]

diffusionCoefficient = 0.2 * 1e-4 * 1e-2; % [ dm^2 s^(-1) ]
molarVolume = 64.12 * 1e-3; % [ dm^3 mol^(-1) ]
%rateCoefficientHydrogen = 0.89 * 1e-2; % [ mol dm^(-2) s^(-1) ]
rateCoefficientTST = 4.5e-4 * 1e-2; % [ mol dm^(-2) s^(-1) ]
equilibriumRateConstant = 10^(-16.5);
inletVelocity = 0.01 * 1e-1; % [ dm s^(-1) ]; Pe = 500
%inletVelocity = 0.001 * 1e-1; % [ dm s^(-1) ]; Pe = 50

initialHydrogenConcentration = 1e-5; % [ mol dm^(-3) ]

dissolutionReactionRate = @(cHydrogen, cCalcium, cCarbonate, cMagnesium) ...
    (cHydrogen^0.5 * rateCoefficientTST) ...
    * (1 - cCalcium * cCarbonate^2 * cMagnesium / equilibriumRateConstant);

% Simulation and discretization parameters (comment for convergence tests)

numberOfSlices = 27; %50;
disp(['numSlices = ', num2str(numberOfSlices)]);
pecletNumber = (inletVelocity * lengthXAxis) / diffusionCoefficient;
fprintf(['Peclet number: ', num2str(pecletNumber), '\n']);

% dt_macro = 1 * 3600; % Pe = 5000
% initialMacroscaleTimeStepSize = 160; % Pe = 50
initialMacroscaleTimeStepSize = 60;
% dt_micro = dt_macro / 100; % 0;
%numMicroscaleTimeSteps = 10;
endTime = 600 * 3600; % [h] Pe = 5000
% endTime = 10 * 60 * 60; % Pe = 50
% endTime = 200 * dt_macro;

%% Computation of time steps in simulation
% stepperType = 'linear';
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
        %         timeStepSizeFactor = 1.5;
        timeStepSizeFactor = 2;
        numTimeSteps = floor(log(endTime / initialMacroscaleTimeStepSize) ...
            /log(timeStepSizeFactor));
        timeSteps = [0; timeStepSizeFactor.^(0:numTimeSteps)'] ...
            * initialMacroscaleTimeStepSize;
        if (timeSteps(end) < endTime)
            timeSteps(end+1) = endTime;
        end
    case 'expmax'
        maximalStep = 0.5 * 10^5;
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
normalSpeedMax = zeros(numberOfSlices, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization of level set method variables

microscaleGrid = FoldedCartesianGrid(dimension, ...
    kron(ones(1, dimension), [-0.5, 0.5]), ...
    numPartitionsMicroscale*ones(1, dimension));
coord = microscaleGrid.coordinates;
coordCell = mat2cell(coord, ones(1, microscaleGrid.nodes), dimension);
helpGridHyPHM = Grid(coord, microscaleGrid.triangles); %same Grid unfolded in HyPHM format

clear coord;

initialLevelSetFunc = @(x) 0.2839 - norm(x);
initialLevelSetDataCells = cell(numberOfSlices, 1);
initialLevelSetDataCells{1} = cellfun(initialLevelSetFunc, coordCell);
for i = 2:numberOfSlices
    initialLevelSetDataCells{i} = initialLevelSetDataCells{1};
end


interfaceNormalVelocity = @(t, x, cHydrogen, cCalcium, cCarbonate, cMagnesium) ...
    spaceScaleFactor * molarVolume * dissolutionReactionRate(cHydrogen, ...
    cCalcium, cCarbonate, cMagnesium);

currentTime = 0;

currentLevelSetDataCells = cell(numberOfSlices, 1);
[currentLevelSetDataCells{:}] = deal(initialLevelSetDataCells{1});
oldLevelSetDataCells = currentLevelSetDataCells;
% TODO Viel zu viele Spalten in levelSetData
% prev: ( ceil( endTime / dt_macro ) + 1 ) columns
levelSetData = NaN(numel(initialLevelSetDataCells{1}), 1);
levelSetData(:, 1) = currentLevelSetDataCells{1};
levelSet = cell(numberOfSlices, 1);
[levelSet{:}] = deal(levelSetData);

clear levelSetData;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization of HyPHM variables

gridHyPHM = domainRectangle(0, lengthXAxis, 0, lengthYAxis, macroscaleStepSize);
macroCoordCell = mat2cell(gridHyPHM.baryT, ones(gridHyPHM.numT, 1), 2);

% timeSteps(end) = []; % At the last time step, no computations are neccessary.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculation of initial effective parameters


diffusionTensors = cell(numberOfSlices, 1);
[diffusionTensors{:}] = deal(NaN(4, numTimeSlices));
permeabilityTensors = cell(numberOfSlices, 1);
[permeabilityTensors{:}] = deal(NaN(4, numTimeSlices));
porosities = cell(numberOfSlices, 1);
[porosities{1:end}] = deal(NaN(numTimeSlices, 1));
clear numTimeSlices;

surfaceArea = cell(numberOfSlices, 1);
%triangleVolumesOld = cell( numberOfSlices, 1 );
timerli = 0;
lil = tic;

for (i = 1:numberOfSlices) %, numProc)
    tic
    [isBalli, radiusEstimate] = isBall(levelSet{i}(:, 1), 0.01); toc
    if false %isBalli
        [D, P, surfaceArea{i}(1), porosities{i}(1)] = DataCircle(radiusEstimate);
        diffusionTensors{i}(:, 1) = D(:);
        permeabilityTensors{i}(:, 1) = P(:) * (meanGrainSize^2);
        %triangleVolumesOld{i}(:,1) = zeros(helpGridHyPHM.numT,1);
    else
        tic
        [cellProblemSystemMatrix, rhs, isDoF, triangleVolumes, triangleSurfaces] ...
            = assembleCellProblem(microscaleGrid, levelSet{i}(:, 1));

        %triangleVolumesOld{i} = triangleVolumes;                                           %reference for adaptive treatment of cell problems
        SOL = solveSystemFE(microscaleGrid, cellProblemSystemMatrix, rhs, isDoF);
        [diffusion, porosities{i}(1)] ...
            = computeDiffusionTensor(microscaleGrid, SOL, triangleVolumes);
        diffusionTensors{i}(:, 1) = diffusion(:);
        timeDiff = timeDiff + toc;
        % porosities{i}(1) = porosity;
        tic
        permeability = computePermeabilityTensor(helpGridHyPHM, levelSet{i}(:, 1));
        permeabilityTensors{i}(:, 1) = permeability(:) * (meanGrainSize^2);
        timePerm = timePerm + toc;
        surfaceArea{i}(1) = sum(triangleSurfaces);
    end

end
timerli = timerli + toc(lil);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization of HyPHM variables (cont.)

flowStepper = Stepper(0:1);
transportStepper = Stepper(timeSteps);

porosityFunc = @(t, x) porosityHelperFun(t, x, porosities, numberOfSlices, 1, lengthXAxis);
diffusionFunc = @(t, x) diffusionHelperFun(t, x, diffusionTensors, ...
    numberOfSlices, 1, lengthXAxis) * diffusionCoefficient;
permeabilityFunc = @(t, x) permeabilityHelperFun(t, x, permeabilityTensors, ...
    numberOfSlices, 1, lengthXAxis);

% Boundary IDs: 1 = down, 2 = right, 3 = up, 4 = left

%Calculate Pressure distribution and Darcy velocity
disp(' ');
disp('Initializing darcy problem...');
phead = Variable(gridHyPHM, flowStepper, 'pressure head', 'P0');
darcyVelocity = Variable(gridHyPHM, flowStepper, 'darcy Velocity', 'RT0');
phead.setdata(0, @(t, x) 0.0);
darcyVelocity.setdata(0, @(t, x) 0.0);

darcy = Transport(gridHyPHM, flowStepper, 'Darcy problem');
darcy.id2D = {2};
darcy.id2F = {1, 3, 4};
darcy.uD.setdata(@(t, x) 0);
darcy.gF.setdata(@(t, x) -inletVelocity*(x(1) < EPS));
darcy.D = Variable(gridHyPHM, flowStepper, 'Permeability', 'P0P0P0P0');
darcy.D.setdata(@(t, x) permeabilityFunc(t, x)./dynamicViscosity);
darcy.F.setdata(@(t, x) 0);
darcy.Q = darcyVelocity;
darcy.U = phead;


flowStepper.next;
tic
darcy.computeLevel('s');
timeDracy = timeDracy + toc;

temp1 = darcy.Q.getdata(1);
flow = Variable(gridHyPHM, flowStepper, 'Flow', 'RT0');
flow.setdata(0, temp1); %variable containing actual flow velocity
%darcy.Q.visualize();
flowStepper = Stepper(0:1);

disp(' ');
disp('Initializing calcium transport...');

calciumConcentration = Variable(gridHyPHM, transportStepper, ...
    'Ca^(2+)', 'P0');
calciumConcentration.setdata(0, @(t, x) 0);
% calciumConcentration.setdata( 0, @(t,x) sqrt( equilibriumRateConstant ) );
calciumTransport = Transport(gridHyPHM, transportStepper, 'Ca^(2+) Transport');
% calciumTransport.id2D = {2};
calciumTransport.id2N = {2};
calciumTransport.id2F = {1, 3, 4};
calciumTransport.U = calciumConcentration;
% calciumTransport.uD = Variable( g_HyPHM, st, 'uD', 'P0E' );
% calciumTransport.uD.setdata( zeros( g_HyPHM.numE, 1 ) );
calciumTransport.gF = Variable(gridHyPHM, transportStepper, 'gF', 'P0E');
calciumTransport.gF.setdata(zeros(gridHyPHM.numE, 1));
calciumTransport.A = Variable(gridHyPHM, transportStepper, 'Porosity', 'P0');
calciumTransport.A.setdata(0, porosityFunc);
calciumTransport.D = Variable(gridHyPHM, transportStepper, 'Diffusion', 'P0P0P0P0');
calciumTransport.D.setdata(0, diffusionFunc);
calciumTransport.C.setdata(0, darcy.Q.getdata(1));
% calciumTransport.C.setdata( P2P2.P2P2toRT0slice( g_HyPHM, flow.getdata(1) ) );
calciumTransport.isUpwind = 'exp';

disp(' ');
disp('Initializing carbonate transport...');

carbonateConcentration = Variable(gridHyPHM, transportStepper, ...
    'CO_3^(2-)', 'P0');
carbonateConcentration.setdata(0, @(t, x) 0);
% carbonateConcentration.setdata( 0, @(t,x) sqrt( equilibriumRateConstant ) );
carbonateTransport = Transport(gridHyPHM, transportStepper, 'CO_3^(2-) Transport');
% carbonateTransport.id2D = {2};
carbonateTransport.id2N = {2};
carbonateTransport.id2F = {1, 3, 4};
carbonateTransport.U = carbonateConcentration;
% carbonateTransport.uD = Variable( g_HyPHM, st, 'uD', 'P0E' );
% carbonateTransport.uD.setdata( zeros( g_HyPHM.numE, 1 ) );
carbonateTransport.gF = Variable(gridHyPHM, transportStepper, 'gF', 'P0E');
carbonateTransport.gF.setdata(zeros(gridHyPHM.numE, 1));
carbonateTransport.A = Variable(gridHyPHM, transportStepper, 'Porosity', 'P0');
carbonateTransport.A.setdata(0, porosityFunc);
carbonateTransport.D = Variable(gridHyPHM, transportStepper, 'Diffusion', 'P0P0P0P0');
carbonateTransport.D.setdata(0, diffusionFunc);
carbonateTransport.C.setdata(0, darcy.Q.getdata(1));
% carbonateTransport.C.setdata( P2P2.P2P2toRT0slice( g_HyPHM, flow.getdata(1) ) );
carbonateTransport.isUpwind = 'exp';

disp(' ');
disp('Initializing magnesium transport...');

magnesiumConcentration = Variable(gridHyPHM, transportStepper, ...
    'Mg^(2+)', 'P0');
magnesiumConcentration.setdata(0, @(t, x) 0);
% magnesiumConcentration.setdata( 0, @(t,x) sqrt( equilibriumRateConstant ) );
magnesiumTransport = Transport(gridHyPHM, transportStepper, 'MG^(2+) Transport');
% magnesiumTransport.id2D = {2};
magnesiumTransport.id2N = {2};
magnesiumTransport.id2F = {1, 3, 4};
magnesiumTransport.U = magnesiumConcentration;
% magnesiumTransport.uD = Variable( g_HyPHM, st, 'uD', 'P0E' );
% magnesiumTransport.uD.setdata( zeros( g_HyPHM.numE, 1 ) );
magnesiumTransport.gF = Variable(gridHyPHM, transportStepper, 'gF', 'P0E');
magnesiumTransport.gF.setdata(zeros(gridHyPHM.numE, 1));
magnesiumTransport.A = Variable(gridHyPHM, transportStepper, 'Porosity', 'P0');
magnesiumTransport.A.setdata(0, porosityFunc);
magnesiumTransport.D = Variable(gridHyPHM, transportStepper, 'Diffusion', 'P0P0P0P0');
magnesiumTransport.D.setdata(0, diffusionFunc);
magnesiumTransport.C.setdata(0, darcy.Q.getdata(1));
% magnesiumTransport.C.setdata( P2P2.P2P2toRT0slice( g_HyPHM, flow.getdata(1) ) );
magnesiumTransport.isUpwind = 'exp';

disp(' ');
disp('Initializing hydrogen transport...');

hydrogenConcentration = Variable(gridHyPHM, transportStepper, ...
    'H^+_pH', 'P0');
hydrogenConcentration.setdata(0, @(t, x) initialHydrogenConcentration);
hydrogenTransport = Transport(gridHyPHM, transportStepper, 'H^+ Transport');
hydrogenTransport.id2N = {2};
hydrogenTransport.id2F = {1, 3, 4};
% hydrogenTransport.id2D = {4};
hydrogenTransport.U = hydrogenConcentration;
hydrogenTransport.gF.setdata( ...
    -initialHydrogenConcentration*inletVelocity*(gridHyPHM.baryE(:, 1) < EPS));
% hydrogenTransport.gF.setdata( @(t,x) -1e-7 );
% hydrogenTransport.uD.setdata( @(t,x) initialHydrogenConcentration );
hydrogenTransport.A.setdata(0, porosityFunc);
hydrogenTransport.C.setdata(0, darcy.Q.getdata(1));
hydrogenTransport.D = Variable(gridHyPHM, transportStepper, 'Diffusion', 'P0P0P0P0');
hydrogenTransport.D.setdata(0, diffusionFunc);
hydrogenTransport.isUpwind = 'exp';

speciesCells = {calciumTransport; carbonateTransport; hydrogenTransport; magnesiumTransport};

preprocessingTime = toc; % Preprocessing

disp(' ');
disp(['Preprocessing done in ', num2str(preprocessingTime), ' seconds.']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    % Level set evolution step

    disp('Evolution of level set ...');
    tic;
    coordCell = mat2cell(microscaleGrid.coordinates, ...
        ones(1, microscaleGrid.nodes), dimension);
    timeCell = num2cell(currentTime*ones(microscaleGrid.nodes, 1));

    calciumSliceAverageCells = cell(numberOfSlices, 1);
    carbonateSliceAverageCells = cell(numberOfSlices, 1);
    magnesiumSliceAverageCells = cell(numberOfSlices, 1);
    hydrogenSliceAverageCells = cell(numberOfSlices, 1);
    flowSliceAverageCells = cell(numberOfSlices, 1);
    calciumData = calciumConcentration.getdata(timeIterationStep-1);
    carbonateData = carbonateConcentration.getdata(timeIterationStep-1);
    magnesiumData = magnesiumConcentration.getdata(timeIterationStep-1);
    hydrogenData = hydrogenConcentration.getdata(timeIterationStep-1);


    data = flow.getdata(0);
    flowVelocitiesCartesian = RT0.RT0toP0P0slice(gridHyPHM, data);

    %average concentration over macro subdomains related to unit cells
    for i = 1:numberOfSlices

        isInArea = ...
            (gridHyPHM.baryT(:, 1) >= ((i - 1) / numberOfSlices) * lengthXAxis) & ...
            (gridHyPHM.baryT(:, 1) < (i / numberOfSlices) * lengthXAxis) & ...
            (gridHyPHM.baryT(:, 2) > lengthYAxis / 26);
        calciumSliceAverageCells{i} = sum(calciumData(isInArea).* ...
            gridHyPHM.areaT(isInArea)) / sum(gridHyPHM.areaT(isInArea));
        carbonateSliceAverageCells{i} = sum(carbonateData(isInArea).* ...
            gridHyPHM.areaT(isInArea)) / sum(gridHyPHM.areaT(isInArea));
        magnesiumSliceAverageCells{i} = sum(magnesiumData(isInArea).* ...
            gridHyPHM.areaT(isInArea)) / sum(gridHyPHM.areaT(isInArea));
        hydrogenSliceAverageCells{i} = sum(hydrogenData(isInArea).* ...
            gridHyPHM.areaT(isInArea)) / sum(gridHyPHM.areaT(isInArea));
        flowSliceAverageCells{i} = sum(flowVelocitiesCartesian(isInArea, :).* ...
            gridHyPHM.areaT(isInArea), 1) / sum(gridHyPHM.areaT(isInArea));
    end


    tic
    for i = 1:numberOfSlices
        %         normalSpeed = cellfun( interfaceNormalVelocity, timeCell, coordCell, ...
        %             num2cell( calciumSliceAverageCells{i} * ones( microscaleGrid.nodes, 1 ) ), ...
        %             num2cell( carbonateSliceAverageCells{i} * ones( microscaleGrid.nodes, 1 ) ) );
        normalSpeed = cellfun(interfaceNormalVelocity, timeCell, coordCell, ...
            num2cell(hydrogenSliceAverageCells{i} ...
            * ones(microscaleGrid.nodes, 1)), ...
            num2cell(calciumSliceAverageCells{i} ...
            * ones(microscaleGrid.nodes, 1)), ...
            num2cell(carbonateSliceAverageCells{i} ...
            * ones(microscaleGrid.nodes, 1)), ...
            num2cell(magnesiumSliceAverageCells{i} ...
            * ones(microscaleGrid.nodes, 1)));

        levelSetPartialDerivatives = zeros(numel(currentLevelSetDataCells{i}), 1);
        for d = 1:dimension
            levelSetPartialDerivatives(:, d) = levelSetPartialDerivative( ...
                microscaleGrid, ...
                currentLevelSetDataCells{i}, d);
        end

        normalSpeedModifier = levelSetPartialDerivatives * flowSliceAverageCells{i}';
        %
        normalSpeedModifier = arrayfun(@normalVelocityModifierFunc, ...
            normalSpeedModifier);
        normalSpeed = normalSpeed .* normalSpeedModifier;
        normalSpeedMax(i) = max(abs(normalSpeed));

        CFL = 1 / 4 * 1 / numPartitionsMicroscale / normalSpeedMax(i);
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

        levelSet{i}(:, timeIterationStep + 1) = currentLevelSetDataCells{i};
    end
    timeLevelSet = timeLevelSet + toc;
    %     t_old = currentTime;
    %     currentTime = currentTime + dt_macro;
    %     oldLevelSetDataCells = currentLevelSetDataCells;
    currentTime = currentTime + macroscaleTimeStepSize;

    levelSetEvolutionTime(timeIterationStep) = toc;
    disp(['    ... done in ', ...
        num2str(levelSetEvolutionTime(timeIterationStep)), ' seconds.']);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Calculation of effective parameters

    disp('Calculation of effective parameters ...');

    count = 0;
    lil = tic;
    %Lev = parallel.pool.Constant(levelSet);
    for (i = 1:numberOfSlices) %, numProc)
        [isBalli, radiusEstimate] = isBall(levelSet{i}(:, timeIterationStep + 1), 0.01);
        if false %isBalli
            [D, P, surfaceArea{i}(timeIterationStep + 1), porosities{i}(timeIterationStep + 1)] = DataCircle(radiusEstimate);
            diffusionTensors{i}(:, timeIterationStep + 1) = D(:);
            permeabilityTensors{i}(:, timeIterationStep + 1) = P(:) * (meanGrainSize^2);
            %triangleVolumesOld{i}(:,timeIterationStep + 1) = zeros(helpGridHyPHM.numT,1);
        else
            tic
            [cellProblemSystemMatrix, rhs, isDoF, triangleVolumes, ...
                triangleSurfaces] = assembleCellProblem(microscaleGrid, ...
                levelSet{i}(:, timeIterationStep + 1));
            timeDiff = timeDiff + toc;
            %if sum(abs(triangleVolumesOld{i}-triangleVolumes)>0.05 * 2 / numPartitionsMicroscale^2)>eps                 % if recalculation necessary
            %    triangleVolumesOld{i} = triangleVolumes;
            if transportStepper.curtau * normalSpeedMax(i) > -0.01 * lengthXAxis
                tic
                SOL = solveSystemFE(microscaleGrid, cellProblemSystemMatrix, rhs, isDoF);
                [diffusion, porosities{i}(timeIterationStep + 1)] ...
                    = computeDiffusionTensor(microscaleGrid, SOL, ...
                    triangleVolumes);
                diffusionTensors{i}(:, timeIterationStep + 1) = diffusion(:);
                %             porosities{i}( timeIterationStep ) = porosity;
                timeDiff = timeDiff + toc;

                tic
                permeability = computePermeabilityTensor(helpGridHyPHM, ...
                    levelSet{i}(:, timeIterationStep + 1));

                permeabilityTensors{i}(:, timeIterationStep + 1) = permeability(:) * (meanGrainSize^2);
                timePerm = timePerm + toc;
            else % copy old data
                diffusionTensors{i}(:, timeIterationStep + 1) = diffusionTensors{i}(:, timeIterationStep);
                porosities{i}(timeIterationStep + 1) = porosities{i}(timeIterationStep);
                permeabilityTensors{i}(:, timeIterationStep + 1) = permeabilityTensors{i}(:, timeIterationStep);
                count = count + 1;
            end
            surfaceArea{i}(timeIterationStep + 1) = sum(triangleSurfaces);
        end

    end
    timerli = timerli + toc(lil);

    disp([num2str(count), ' recalculation(s) saved by adaptivity :-)']);
    cellProblemTime(timeIterationStep) = toc;
    disp(['    ... done in ', num2str(cellProblemTime(timeIterationStep)), ...
        ' seconds.']);

    % calculate Darcy velocity and pressure field
    porosityFunc = @(t, x) porosityHelperFun(t, x, porosities, numberOfSlices, timeIterationStep+1, lengthXAxis);
    permeabilityFunc = @(t, x) permeabilityHelperFun(t, x, permeabilityTensors, ...
        numberOfSlices, timeIterationStep+1, lengthXAxis);

    tic
    disp(['Solving Darcy equation...']);
    darcy.D.setdata(1, @(t, x) permeabilityFunc(t, x)./dynamicViscosity);
    flowStepper.next;
    darcy.computeLevel('s');
    timeDracy = timeDracy + toc;
    disp(['Preprocessing done in ', num2str(toc), ' seconds.']);

    temp1 = darcy.Q.getdata(1);

    flow.setdata(0, temp1);
    flowStepper = Stepper(0:1);

    % set advective velocity field
    magnesiumTransport.C.setdata(timeIterationStep, darcy.Q.getdata(1));
    hydrogenTransport.C.setdata(timeIterationStep, darcy.Q.getdata(1));
    calciumTransport.C.setdata(timeIterationStep, darcy.Q.getdata(1));
    carbonateTransport.C.setdata(timeIterationStep, darcy.Q.getdata(1));


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Macroscopic transport step

    porosityFunc = @(t, x) porosityHelperFun(t, x, porosities, numberOfSlices, timeIterationStep+1, lengthXAxis);
    diffusionFunc = @(t, x) diffusionHelperFun(t, x, diffusionTensors, ...
        numberOfSlices, timeIterationStep+1, lengthXAxis) * diffusionCoefficient;
    permeabilityFunc = @(t, x) permeabilityHelperFun(t, x, permeabilityTensors, ...
        numberOfSlices, timeIterationStep+1, lengthXAxis);


    %     surfaceAreaFun = @(x) areaHelperFun( -1, x, microscaleGrid, levelSet, numberOfSlices, timeStep );
    surfaceAreaFunc = @(x) areaHelperFun(-1, x, surfaceArea, numberOfSlices, ...
        timeIterationStep+1, lengthXAxis);


    specificSurfaceArea = surfaceAreaFunc(gridHyPHM.baryT') * spaceScaleFactor;
    % TODO evaluate dissolutionReactionRate and other related function
    % handles (really needed?)

    % nonlinear reaction rates
    nonlinearFunc = cell(4, 1);
    nonlinearFunc{1} = @(x, y, z, w) specificSurfaceArea ...
        .* (sqrt(z) * rateCoefficientTST) ...
        .* (1 - x .* y .* y .* w ./ equilibriumRateConstant);
    nonlinearFunc{2} = @(x, y, z, w) 2 * nonlinearFunc{1}(x, y, z, w);
    nonlinearFunc{3} = @(x, y, z, w) 0;
    nonlinearFunc{4} = @(x, y, z, w) nonlinearFunc{1}(x, y, z, w);

    nonlinearJacFunc = cell(4, 4);
    nonlinearJacFunc{1, 1} = @(x, y, z, w) specificSurfaceArea ...
        .* (sqrt(z) * rateCoefficientTST) ...
        .* (-y .* y .* w ./ equilibriumRateConstant);
    nonlinearJacFunc{1, 2} = @(x, y, z, w) specificSurfaceArea ...
        .* (sqrt(z) * rateCoefficientTST) ...
        .* (-2 * x .* y .* w ./ equilibriumRateConstant);
    nonlinearJacFunc{1, 3} = @(x, y, z, w) specificSurfaceArea ...
        .* (0.5 * 1 ./ sqrt(z) * rateCoefficientTST) ...
        .* (1 - x .* y .* y .* w ./ equilibriumRateConstant);
    nonlinearJacFunc{1, 4} = @(x, y, z, w) specificSurfaceArea ...
        .* (sqrt(z) * rateCoefficientTST) ...
        .* (-x .* y .* y ./ equilibriumRateConstant);

    nonlinearJacFunc{2, 1} = @(x, y, z, w) 2 * nonlinearJacFunc{1, 1}(x, y, z, w);
    nonlinearJacFunc{2, 2} = @(x, y, z, w) 2 * nonlinearJacFunc{1, 2}(x, y, z, w);
    nonlinearJacFunc{2, 3} = @(x, y, z, w) 2 * nonlinearJacFunc{1, 3}(x, y, z, w);
    nonlinearJacFunc{2, 4} = @(x, y, z, w) 2 * nonlinearJacFunc{1, 4}(x, y, z, w);

    nonlinearJacFunc{3, 1} = @(x, y, z, w) 0;
    nonlinearJacFunc{3, 2} = @(x, y, z, w) 0;
    nonlinearJacFunc{3, 3} = @(x, y, z, w) 0;
    nonlinearJacFunc{3, 4} = @(x, y, z, w) 0;

    nonlinearJacFunc{4, 1} = @(x, y, z, w) nonlinearJacFunc{1, 1}(x, y, z, w);
    nonlinearJacFunc{4, 2} = @(x, y, z, w) nonlinearJacFunc{1, 2}(x, y, z, w);
    nonlinearJacFunc{4, 3} = @(x, y, z, w) nonlinearJacFunc{1, 3}(x, y, z, w);
    nonlinearJacFunc{4, 4} = @(x, y, z, w) nonlinearJacFunc{1, 4}(x, y, z, w);

    % set parameters for transport instances
    calciumTransport.A.setdata(timeIterationStep, porosityFunc);
    calciumTransport.D.setdata(timeIterationStep, diffusionFunc);
    carbonateTransport.A.setdata(timeIterationStep, porosityFunc);
    carbonateTransport.D.setdata(timeIterationStep, diffusionFunc);
    hydrogenTransport.A.setdata(timeIterationStep, porosityFunc);
    hydrogenTransport.D.setdata(timeIterationStep, diffusionFunc);
    magnesiumTransport.A.setdata(timeIterationStep, porosityFunc);
    magnesiumTransport.D.setdata(timeIterationStep, diffusionFunc);

    tic
    newtonIteration(speciesCells, nonlinearFunc, nonlinearJacFunc, 10);

    timeNewton = timeNewton + toc;

end % while
toc(all)
timeNewton
timeLevelSet
timeDiff
timePerm
timeDracy
PermeabilityTensorTest = permeabilityTensors{1};
CalciumConcentrationTest = calciumTransport.U.getdata(17);
GridTest = gridHyPHM;


% calciumConcentration.visualize;
% carbonateConcentration.visualize;
% hydrogenConcentration.visualize;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Helper functions

function por = porosityHelperFun(t, x, porosityCells, numberOfSlices, timeStep, lengthXAxis)
%return porosity value for each macroscopic triangle
por = 0;
x = x ./ lengthXAxis;
for i = 1:numberOfSlices - 1
    por = por + porosityCells{i}(timeStep) * ...
        (x(1) >= (i - 1) / numberOfSlices & x(1) < i / numberOfSlices);
end
por = por + porosityCells{numberOfSlices}(timeStep) * ...
    (x(1) >= (numberOfSlices - 1) / numberOfSlices & x(1) <= 1);
%     por = porosityCells{ findSlice( x, numberOfSlices ) }( timeStep );

end


function perm = permeabilityHelperFun(t, x, permeabilityCells, numberOfSlices, timeStep, lengthXAxis)
%return permeability tensor for each macroscopic triangle
perm = zeros(size(x, 2), 4);
x = x ./ lengthXAxis;
for i = 1:numberOfSlices - 1
    perm = perm + permeabilityCells{i}(:, timeStep)' .* ...
        (x(1, :) >= (i - 1) / numberOfSlices & x(1, :) < i / numberOfSlices)';
end
perm = perm + permeabilityCells{numberOfSlices}(:, timeStep)' .* ...
    (x(1, :) >= (numberOfSlices - 1) / numberOfSlices & x(1, :) <= 1)';
end

function diff = diffusionHelperFun(t, x, diffusionCells, numberOfSlices, timeStep, lengthXAxis)
%return diffusion tensor for each macroscopic triangle
diff = zeros(size(x, 2), 4);
x = x ./ lengthXAxis;
for i = 1:numberOfSlices - 1
    diff = diff + diffusionCells{i}(:, timeStep)' .* ...
        (x(1, :) >= (i - 1) / numberOfSlices & x(1, :) < i / numberOfSlices)';
end

diff = diff + diffusionCells{numberOfSlices}(:, timeStep)' .* ...
    (x(1, :) >= (numberOfSlices - 1) / numberOfSlices & x(1, :) <= 1)';
%     diff = reshape( diffusionCells{ findSlice( x, numberOfSlices ) }( :, ...
%         timeStep ), 2, 2 );

end

function area = areaHelperFun(t, x, surfaceAreaCells, numberOfSlices, timeStep, lengthXAxis)
%return surface area value for each macroscopic triangle
area = zeros(size(x, 2), 1);
x(1, :) = x(1, :) ./ lengthXAxis;
for i = 1:numberOfSlices
    area = area + surfaceAreaCells{i}(timeStep) * ...
        (x(1, :) >= (i - 1) / numberOfSlices & x(1, :) < i / numberOfSlices)';
end

end

% function area = areaHelperFun( t, x, grid, levelSetCells, numberOfSlices, timeStep )
%
% %     area = 0;
% %     for i = 1:numberOfSlices
% %         area = area + surfaceIntegral( grid, levelSetCells{i}( :, timeStep ), ...
% %             @(x) 1 ) * ...
% %             ( x(1) >= (i-1)/numberOfSlices & x(1) < i/numberOfSlices );
% %     end
%
%     area = surfaceIntegral( grid, ...
%         levelSetCells{ findSlice( x, numberOfSlices ) }( :, timeStep ), @(x) 1 );
%
% end


function sliceNumber = findSlice(x, numberOfSlices, lengthXAxis)
% identify number of slice to which coordinates x belong
sliceNumber = 1 + fix(x(1)./lengthXAxis*numberOfSlices);
if sliceNumber > numberOfSlices
    sliceNumber = numberOfSlices;
end
end


function y = normalVelocityModifierFunc(x)
%modify uniform normal velocity in order to obtain droplet shapes
%grains
positiveLimitValue = 0.75;
%expFactor = 5e+3; % Value for convergence tests
expFactor = 5e+3;
positive = (x > 0.0);
y = exp(-expFactor*x.^2);
y(positive) = y(positive) + positiveLimitValue ...
    * (1.0 - exp(-expFactor * x(positive).^2));

end