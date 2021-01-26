%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% General setting variables
global Solver
Solver = 'BlockPrec'
all = tic;
global EPS;
EPS = eps;
global timeNewton;
global timeLevelSet;
global timeDiff;
global timePerm;
global timeDarcy;
global numPartitionsMicroscale;
global macroscaleStepSize;
global PermeabilityTensorTest;
global CalciumConcentrationTest;
global GridTest;

numProc = 1;
lengthXAxis = 1 * 0.1; % [ dm ]
lengthYAxis = 0.5 * 0.1; % [ dm ]
numPartitionsMicroscale = 100; % Number of partitions in each direction
numPartitionsMicroscaleCoarse = 50;
macroscaleStepSize = lengthXAxis / 100; %(27*4); % [ dm ]
meanGrainSize = 10^(-4); % [ dm ]


tic; % Preprocessing

% Physical parameters
dimension = 2;
spaceScaleFactor = 10 * 10; %10*sqrt(13*27); % spaceScaleFactor [length of Y] = 1 [cm]

% Parameters in dm and s
dynamicViscosity = 1100e-7; % [ kg dm^(-1) s^(-1) ]

diffusionCoefficient = 0.2 * 1e-4 * 1e-2; % [ dm^2 s^(-1) ]

molarVolume = cell(2, 1);
molarVolume{1} = 39.63 * 1e-3;
molarVolume{2} = 64.12 * 1e-3; % [ dm^3 mol^(-1) ]

rateCoefficientHydrogen = 0.89 * 1e-2; % [ mol dm^(-2) s^(-1) ]

rateCoefficientTSTC = 6.6e-7 * 1e-2; % [ mol dm^(-2) s^(-1) ]
rateCoefficientTSTD = 4.5e-4 * 1e-2; % [ mol dm^(-2) s^(-1) ]

equilibriumRateConstant = cell(2, 1);
equilibriumRateConstant{1} = 10^(-8.234);
equilibriumRateConstant{2} = 10^(-16.5);

inletVelocity = 0.02 * 1e-1; % [ dm s^(-1) ]; Pe = 500
%inletVelocity = 0.001 * 1e-1; % [ dm s^(-1) ]; Pe = 50

initialHydrogenConcentration = 1e-5; % [ mol dm^(-3) ]

dissolutionReactionRate = cell(2, 1);
dissolutionReactionRate{1} = @(cHydrogen, cCalcium, cCarbonate, cMagnesium) ...
    (rateCoefficientHydrogen * cHydrogen + rateCoefficientTSTC) ...
    * (1 - cCalcium * cCarbonate / equilibriumRateConstant{1});

dissolutionReactionRate{2} = @(cHydrogen, cCalcium, cCarbonate, cMagnesium) ...
    (cHydrogen^0.5 * rateCoefficientTSTD) ...
    * (1 - cCalcium * cCarbonate^2 * cMagnesium / equilibriumRateConstant{2});

% Simulation and discretization parameters (comment for convergence tests)

numberOfSlices = 10; %50;
microRunSinceUpdate = zeros(numberOfSlices, 1);
disp(['numSlices = ', num2str(numberOfSlices)]);
pecletNumber = (inletVelocity * lengthXAxis) / diffusionCoefficient;
fprintf(['Peclet number: ', num2str(pecletNumber), '\n']);

% dt_macro = 1 * 3600; % Pe = 5000
% initialMacroscaleTimeStepSize = 160; % Pe = 50
initialMacroscaleTimeStepSize = 60 * 2^(-3);
% dt_micro = dt_macro / 100; % 0;
%numMicroscaleTimeSteps = 10;
endTime = 1000 * 3600; % [h] Pe = 5000
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
        maximalStep = 0.2 * 10^5;
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

microscaleGridCoarse = FoldedCartesianGrid(dimension, ...
    kron(ones(1, dimension), [-0.5, 0.5]), ...
    numPartitionsMicroscaleCoarse*ones(1, dimension));
coord = microscaleGridCoarse.coordinates;
coordCell = mat2cell(coord, ones(1, microscaleGridCoarse.nodes), dimension);
helpGridHyPHM = Grid(coord, microscaleGridCoarse.triangles); %same Grid unfolded in HyPHM format

microscaleGrid = FoldedCartesianGrid(dimension, ...
    kron(ones(1, dimension), [-0.5, 0.5]), ...
    numPartitionsMicroscale*ones(1, dimension));
coord = microscaleGrid.coordinates;
coordCell = mat2cell(coord, ones(1, microscaleGrid.nodes), dimension);
clear coord;

microCells = initUnitCells(microscaleGrid, numberOfSlices, numTimeSlices, true);
interfaceNormalVelocity = cell(2, 1);
for j = 1:2
    interfaceNormalVelocity{j} = @(t, x, cHydrogen, cCalcium, cCarbonate, cMagnesium) ...
        spaceScaleFactor * molarVolume{j} * dissolutionReactionRate{j}(cHydrogen, ...
        cCalcium, cCarbonate, cMagnesium);
end

currentTime = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization of HyPHM variables

gridHyPHM = domainRectangle(0, lengthXAxis, 0, lengthYAxis, macroscaleStepSize);
macroCoordCell = mat2cell(gridHyPHM.baryT, ones(gridHyPHM.numT, 1), 2);

% timeSteps(end) = []; % At the last time step, no computations are neccessary.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculation of initial effective parameters

% diffusionTensorsData = NaN( 4, numTimeSlices );
% diffusionTensors = cell( numberOfSlices, 1 );
% [ diffusionTensors{1:end} ] = deal( diffusionTensorsData );
% porosities = cell( numberOfSlices, 1 );
% [ porosities{1:end} ] = deal( porositiesData );
diffusionTensors = cell(numberOfSlices, 1);
[diffusionTensors{:}] = deal(NaN(4, numTimeSlices));
permeabilityTensors = cell(numberOfSlices, 1);
[permeabilityTensors{:}] = deal(NaN(4, numTimeSlices));
porosities = cell(numberOfSlices, 1);
[porosities{1:end}] = deal(NaN(numTimeSlices, 1));

surfaceArea = cell(numberOfSlices, 1);
for j = 1:numberOfSlices
    surfaceArea{j} = cell(2, 1);
    [surfaceArea{j}{:}] = deal(NaN(numTimeSlices, 1));
end
clear numTimeSlices;
%triangleVolumesOld = cell( numberOfSlices, 1 );
timerli = 0;
hydroData = tic;
[Xcoarse, Ycoarse] = meshgrid(linspace(-0.5, 0.5, numPartitionsMicroscaleCoarse + 1));
[X, Y] = meshgrid(linspace(-0.5, 0.5, numPartitionsMicroscale + 1));

for (i = 1:numberOfSlices) %, numProc)
    tic
    %    [isBalli, radiusEstimate] = isBall(levelSet{i}(:,1), 0.01);toc
    if false %isBalli
        %         [D,P,surfaceArea{i}(1), porosities{i}(1)] = DataCircle(radiusEstimate);
        %         diffusionTensors{i}(:,1) = D(:);
        %         permeabilityTensors{i}(:,1) = P(:) * (meanGrainSize ^ 2);
        %         %triangleVolumesOld{i}(:,1) = zeros(helpGridHyPHM.numT,1);
    else
        tic

        coarseLevel = interp2(X, Y, reshape(microCells.signed{i}(:, 1), numPartitionsMicroscale + 1, numPartitionsMicroscale + 1), Xcoarse, Ycoarse);

        [cellProblemSystemMatrix, rhs, isDoF, triangleVolumes, triangleSurfaces] ...
            = assembleCellProblem(microscaleGridCoarse, coarseLevel);

        %triangleVolumesOld{i} = triangleVolumes;                                           %reference for adaptive treatment of cell problems
        SOL = solveSystemFE(microscaleGridCoarse, cellProblemSystemMatrix, rhs, isDoF);
        [diffusion, porosities{i}(1)] ...
            = computeDiffusionTensor(microscaleGridCoarse, SOL, triangleVolumes);
        diffusionTensors{i}(:, 1) = diffusion(:);
        timeDiff = timeDiff + toc;
        % porosities{i}(1) = porosity;
        tic
        permeability = computePermeabilityTensor(helpGridHyPHM, coarseLevel);
        permeabilityTensors{i}(:, 1) = permeability(:) * (meanGrainSize^2);
        timePerm = timePerm + toc;
        surfaceArea{i}{1}(1) = microCells.interfaceLength{i}(1);
        surfaceArea{i}{2}(1, 1) = microCells.interfaceLength{i}(1, 2);
    end

end
timerli = timerli + toc(hydroData);
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
timeDarcy = timeDarcy + toc;

temp1 = darcy.Q.getdata(1);
temp2 = temp1;
for ed = 1:gridHyPHM.numE
    temp2(ed) = temp1(ed); % ./sqrt(1-porosityFunc(0,mean(gridHyPHM.coordV(gridHyPHM.V0E(ed,:),:))));
end
flow = Variable(gridHyPHM, flowStepper, 'Flow', 'RT0');
flow.setdata(0, temp2); %variable containing actual flow velocity
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
%magnesiumTransport.id2D = {4};
magnesiumTransport.id2N = {2};
magnesiumTransport.id2F = {1, 3, 4};
magnesiumTransport.U = magnesiumConcentration;
% magnesiumTransport.uD = Variable( g_HyPHM, st, 'uD', 'P0E' );
% magnesiumTransport.uD.setdata( 10^(-2)*ones( gridHyPHM.numE, 1 ) );
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

% calciumDataOld = calciumConcentration.getdata( 0 );
% carbonateDataOld = carbonateConcentration.getdata( 0 );
% hydrogenDataOld = hydrogenConcentration.getdata( 0 );

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

    for i = 1:numberOfSlices

        isInArea = ...
            (gridHyPHM.baryT(:, 1) >= ((i - 1) / numberOfSlices) * lengthXAxis) & ...
            (gridHyPHM.baryT(:, 1) < (i / numberOfSlices) * lengthXAxis);
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
    newlol = tic;

    normalSpeeds = cell(numberOfSlices, 1);
    normalSpeedMax = zeros(numberOfSlices, 1);

    for i = 1:numberOfSlices

        normalSpeed = cell(2, 1);
        for j = 1:2
            normalSpeed{j} = interfaceNormalVelocity{j}(currentTime, [0, 0], ...
                hydrogenSliceAverageCells{i}, ...
                calciumSliceAverageCells{i}, ...
                carbonateSliceAverageCells{i}, ...
                magnesiumSliceAverageCells{i});
        end

        normalSpeed = [0, sign(microCells.interfaceLength{i}(microCells.memoryStep, 1)) * normalSpeed{1}, sign(microCells.interfaceLength{i}(microCells.memoryStep, 2)) * normalSpeed{2}; ...
            - sign(microCells.interfaceLength{i}(microCells.memoryStep, 1)) * normalSpeed{1}, 0, 0; -sign(microCells.interfaceLength{i}(microCells.memoryStep, 2)) * normalSpeed{2}, 0, 0];
        normalSpeeds{i} = normalSpeed;
        normalSpeedMax(i) = max(abs(normalSpeed(:)));


        %microscaleTimeStepSize = transportStepper.curtau / numMicroscaleTimeSteps;
    end

    microCells = evolveStep(microCells, normalSpeeds, ...
        timeIterationStep, transportStepper.curtau, numPartitionsMicroscale);

    toc(newlol)
    timeLevelSet = timeLevelSet + toc;
    %     t_old = currentTime;
    %     currentTime = currentTime + dt_macro;
    %     oldLevelSetDataCells = currentLevelSetDataCells;
    currentTime = currentTime + macroscaleTimeStepSize;

    levelSetEvolutionTime(timeIterationStep) = toc;
    disp(['    ... done in ', ...
        num2str(levelSetEvolutionTime(timeIterationStep)), ' seconds.']);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Plotting of interface

    %     PHI = g.reshape( phi2 );
    %     set( hand, 'Zdata', PHI, 'LineColor', 'b' );
    %     drawnow;
    % %     saveas( fig, [ 'data/MultipleCellProblems/Triangle', ...
    % %              num2str( timeStep ), '.png' ] );
    %     if ( 60 < timeStep && timeStep <= 70 )
    %         print( fig, '-dpng', [ 'data/MultipleCellProblems/Triangle', ...
    %              num2str( timeStep ), '.png' ] );
    %     end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Calculation of effective parameters

    disp('Calculation of effective parameters ...');

    count = 0;
    Balls = 0;
    hydroData = tic;
    %Lev = parallel.pool.Constant(levelSet);

    for (i = 1:numberOfSlices) %, numProc)
        microRunSinceUpdate(i) = microRunSinceUpdate(i) + transportStepper.curtau * normalSpeedMax(i);
        %         [isBalli, radiusEstimate] = isBall(levelSet{i}(:,timeIterationStep + 1), 0.01);
        if false %isBalli
            %             [D,P,surfaceArea{i}(timeIterationStep + 1), porosities{i}(timeIterationStep + 1)] = DataCircle(radiusEstimate);
            %             diffusionTensors{i}(:,timeIterationStep + 1) = D(:);
            %             permeabilityTensors{i}(:,timeIterationStep + 1) = P(:) * (meanGrainSize ^ 2);
            %             %triangleVolumesOld{i}(:,timeIterationStep + 1) = zeros(helpGridHyPHM.numT,1);
            %             Balls = Balls + 1;
        else
            tic
            coarseLevel = interp2(X, Y, reshape(microCells.signed{i}(:, microCells.memoryStep), numPartitionsMicroscale + 1, numPartitionsMicroscale + 1), Xcoarse, Ycoarse);
            [cellProblemSystemMatrix, rhs, isDoF, triangleVolumes, ...
                triangleSurfaces] = assembleCellProblem(microscaleGridCoarse, ...
                coarseLevel);
            timeDiff = timeDiff + toc;

            if microRunSinceUpdate(i) > 0.1 / numPartitionsMicroscaleCoarse
                tic
                microRunSinceUpdate(i) = 0;
                SOL = solveSystemFE(microscaleGridCoarse, cellProblemSystemMatrix, rhs, isDoF);
                [diffusion, porosities{i}(timeIterationStep + 1)] ...
                    = computeDiffusionTensor(microscaleGridCoarse, SOL, ...
                    triangleVolumes);
                diffusionTensors{i}(:, timeIterationStep + 1) = diffusion(:);
                %             porosities{i}( timeIterationStep ) = porosity;
                timeDiff = timeDiff + toc;

                tic
                permeability = computePermeabilityTensor(helpGridHyPHM, ...
                    coarseLevel);

                permeabilityTensors{i}(:, timeIterationStep + 1) = permeability(:) * (meanGrainSize^2);
                timePerm = timePerm + toc;
            else % copy old data
                diffusionTensors{i}(:, timeIterationStep + 1) = diffusionTensors{i}(:, timeIterationStep);
                porosities{i}(timeIterationStep + 1) = porosities{i}(timeIterationStep);
                permeabilityTensors{i}(:, timeIterationStep + 1) = permeabilityTensors{i}(:, timeIterationStep);
                count = count + 1;
            end
            surfaceArea{i}{1}(timeIterationStep + 1, 1) = microCells.interfaceLength{i}(microCells.memoryStep, 1);
            surfaceArea{i}{2}(timeIterationStep + 1, 1) = microCells.interfaceLength{i}(microCells.memoryStep, 2);
        end

    end

    timerli = timerli + toc(hydroData);

    disp([num2str(Balls), ' recalculation(s) saved since circular solid :-)']);
    disp([num2str(count), ' recalculation(s) saved by adaptivity :-)']);
    cellProblemTime(timeIterationStep) = toc;
    disp(['    ... done in ', num2str(cellProblemTime(timeIterationStep)), ...
        ' seconds.']);


    porosityFunc = @(t, x) porosityHelperFun(t, x, porosities, numberOfSlices, timeIterationStep+1, lengthXAxis);
    permeabilityFunc = @(t, x) permeabilityHelperFun(t, x, permeabilityTensors, ...
        numberOfSlices, timeIterationStep+1, lengthXAxis);

    tic
    disp(['Solving Darcy equation...']);
    darcy.D.setdata(1, @(t, x) permeabilityFunc(t, x)./dynamicViscosity);
    flowStepper.next;
    darcy.computeLevel('s');
    timeDarcy = timeDarcy + toc;
    disp(['Preprocessing done in ', num2str(toc), ' seconds.']);

    temp1 = darcy.Q.getdata(1);
    for ed = 1:gridHyPHM.numE
        temp2(ed) = temp1(ed); %./sqrt(1-porosityFunc(0,mean(gridHyPHM.coordV(gridHyPHM.V0E(ed,:),:))));
    end

    flow.setdata(0, temp2);
    flowStepper = Stepper(0:1);

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

    for j = 1:2
        temp = cell(numberOfSlices, 1);
        for k = 1:numberOfSlices
            temp{k} = surfaceArea{k}{j};
        end
        surfaceAreaFunc{j} = @(x) areaHelperFun(-1, x, temp, numberOfSlices, ...
            timeIterationStep+1, lengthXAxis);
        %surfaceAreaFunc{j} = @(x) areaHelperFun( -1, x, {surfaceArea{:,j}}, numberOfSlices, ...
        %    timeIterationStep + 1, lengthXAxis );
        specificSurfaceArea{j} = cellfun(surfaceAreaFunc{j}, macroCoordCell) * spaceScaleFactor;
        clear temp
    end


    % TODO evaluate dissolutionReactionRate and other related function
    % handles (really needed?)
    %  specificSurfaceArea{2} =0*specificSurfaceArea{2};

    nonlinearFunc = cell(4, 1);
    nonlinearFunc{1} = @(x, y, z, w) 1 ...
        .* (sqrt(max(z, eps)) * rateCoefficientTSTD) ...
        .* (1 - x .* y .* y .* w ./ equilibriumRateConstant{2}) .* specificSurfaceArea{2} ...
        +(rateCoefficientHydrogen * z + rateCoefficientTSTC) .* ...
        (1 - x .* y ./ equilibriumRateConstant{1}) .* specificSurfaceArea{1};

    nonlinearFunc{2} = @(x, y, z, w) 2 ...
        .* (sqrt(max(z, eps)) * rateCoefficientTSTD) ...
        .* (1 - x .* y .* y .* w ./ equilibriumRateConstant{2}) .* specificSurfaceArea{2} ...
        +(rateCoefficientHydrogen * z + rateCoefficientTSTC) .* ...
        (1 - x .* y ./ equilibriumRateConstant{1}) .* specificSurfaceArea{1};
    nonlinearFunc{3} = @(x, y, z, w) -1 ...
        .* rateCoefficientHydrogen .* z ...
        .* (1 - x .* y ./ equilibriumRateConstant{1}) .* specificSurfaceArea{1};
    nonlinearFunc{4} = @(x, y, z, w) 1 ...
        .* (sqrt(max(z, eps)) .* rateCoefficientTSTD) ...
        .* (1 - x .* y .* y .* w ./ equilibriumRateConstant{2}) .* specificSurfaceArea{2};

    nonlinearJacFunc = cell(4, 4);
    nonlinearJacFunc{1, 1} = @(x, y, z, w) 1 ...
        .* (sqrt(max(z, eps)) * rateCoefficientTSTD) ...
        .* (-y .* y .* w ./ equilibriumRateConstant{2}) .* specificSurfaceArea{2} ...
        +(rateCoefficientHydrogen * z + rateCoefficientTSTC) .* ...
        (-y ./ equilibriumRateConstant{1}) .* specificSurfaceArea{1};
    nonlinearJacFunc{1, 2} = @(x, y, z, w) 1 ...
        .* (sqrt(max(z, eps)) * rateCoefficientTSTD) ...
        .* (-2 * x .* y .* w ./ equilibriumRateConstant{2}) .* specificSurfaceArea{2} ...
        +(rateCoefficientHydrogen * z + rateCoefficientTSTC) .* ...
        (-x ./ equilibriumRateConstant{1}) .* specificSurfaceArea{1};
    nonlinearJacFunc{1, 3} = @(x, y, z, w) 1 ...
        .* (0.5 ./ sqrt(max(z, eps)) * rateCoefficientTSTD) ...
        .* (1 - x .* y .* y .* w ./ equilibriumRateConstant{2}) .* specificSurfaceArea{2} ...
        +(rateCoefficientHydrogen) .* ...
        (1 - x .* y ./ equilibriumRateConstant{1}) .* specificSurfaceArea{1};
    nonlinearJacFunc{1, 4} = @(x, y, z, w) 1 ...
        .* (sqrt(max(z, eps)) * rateCoefficientTSTD) ...
        .* (-x .* y .* y ./ equilibriumRateConstant{2}) .* specificSurfaceArea{2};

    nonlinearJacFunc{2, 1} = @(x, y, z, w) 2 ...
        .* (sqrt(max(z, eps)) * rateCoefficientTSTD) ...
        .* (-y .* y .* w ./ equilibriumRateConstant{2}) .* specificSurfaceArea{2} ...
        +(rateCoefficientHydrogen * z + rateCoefficientTSTC) .* ...
        (-y ./ equilibriumRateConstant{1}) .* specificSurfaceArea{1};
    nonlinearJacFunc{2, 2} = @(x, y, z, w) 2 ...
        .* (sqrt(max(z, eps)) * rateCoefficientTSTD) ...
        .* (-2 * x .* y .* w ./ equilibriumRateConstant{2}) .* specificSurfaceArea{2} ...
        +(rateCoefficientHydrogen * z + rateCoefficientTSTC) .* ...
        (-x ./ equilibriumRateConstant{1}) .* specificSurfaceArea{1};
    nonlinearJacFunc{2, 3} = @(x, y, z, w) 2 ...
        .* (0.5 ./ sqrt(max(z, eps)) * rateCoefficientTSTD) ...
        .* (1 - x .* y .* y .* w ./ equilibriumRateConstant{2}) .* specificSurfaceArea{2} ...
        +(rateCoefficientHydrogen) .* ...
        (1 - x .* y ./ equilibriumRateConstant{1}) .* specificSurfaceArea{1};
    nonlinearJacFunc{2, 4} = @(x, y, z, w) 2 ...
        .* (sqrt(max(z, eps)) * rateCoefficientTSTD) ...
        .* (-x .* y .* y ./ equilibriumRateConstant{2}) .* specificSurfaceArea{2};

    nonlinearJacFunc{3, 1} = @(x, y, z, w) -1 ...
        .* rateCoefficientHydrogen .* z ...
        .* (-y ./ equilibriumRateConstant{1}) .* specificSurfaceArea{1};
    nonlinearJacFunc{3, 2} = @(x, y, z, w) -1 ...
        .* rateCoefficientHydrogen .* z ...
        .* (-x ./ equilibriumRateConstant{1}) .* specificSurfaceArea{1};
    nonlinearJacFunc{3, 3} = @(x, y, z, w) -1 ...
        .* rateCoefficientHydrogen ...
        .* (1 - x .* y ./ equilibriumRateConstant{1}) .* specificSurfaceArea{1};
    nonlinearJacFunc{3, 4} = @(x, y, z, w) 0 .* z;

    nonlinearJacFunc{4, 1} = @(x, y, z, w) 1 ...
        .* (sqrt(max(z, eps)) * rateCoefficientTSTD) ...
        .* (-y .* y .* w ./ equilibriumRateConstant{2}) .* specificSurfaceArea{2};
    nonlinearJacFunc{4, 2} = @(x, y, z, w) 1 ...
        .* (sqrt(max(z, eps)) * rateCoefficientTSTD) ...
        .* (-2 * x .* y .* w ./ equilibriumRateConstant{2}) .* specificSurfaceArea{2};
    nonlinearJacFunc{4, 3} = @(x, y, z, w) 1 ...
        .* (0.5 ./ sqrt(max(z, eps)) * rateCoefficientTSTD) ...
        .* (1 - x .* y .* y .* w ./ equilibriumRateConstant{2}) .* specificSurfaceArea{2};
    nonlinearJacFunc{4, 4} = @(x, y, z, w) 1 ...
        .* (sqrt(max(z, eps)) * rateCoefficientTSTD) ...
        .* (-x .* y .* y ./ equilibriumRateConstant{2}) .* specificSurfaceArea{2};

    calciumTransport.A.setdata(timeIterationStep, porosityFunc);
    calciumTransport.D.setdata(timeIterationStep, diffusionFunc);
    %calciumTransport.F.setdata( timeIterationStep, calciumTransportRhsData );
    carbonateTransport.A.setdata(timeIterationStep, porosityFunc);
    carbonateTransport.D.setdata(timeIterationStep, diffusionFunc);
    %carbonateTransport.F.setdata( timeIterationStep, carbonateTransportRhsData );
    hydrogenTransport.A.setdata(timeIterationStep, porosityFunc);
    hydrogenTransport.D.setdata(timeIterationStep, diffusionFunc);
    %hydrogenTransport.F.setdata( timeIterationStep, hydrogenTransportRhsData );
    magnesiumTransport.A.setdata(timeIterationStep, porosityFunc);
    magnesiumTransport.D.setdata(timeIterationStep, diffusionFunc);
    %magnesiumTransport.F.setdata( timeIterationStep, magnesiumTransportRhsData );

    tic
    newtonIteration(speciesCells, nonlinearFunc, nonlinearJacFunc, 10);

    timeNewton = timeNewton + toc;

end % while
toc(all)
timeNewton
timeLevelSet
timeDiff
timePerm
timeDarcy
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
por = zeros(size(x, 2), 1);
x = x ./ lengthXAxis;
for i = 1:numberOfSlices - 1
    por = por + porosityCells{i}(timeStep) * ...
        (x(1, :) >= (i - 1) / numberOfSlices & x(1, :) < i / numberOfSlices)';
end
por = por + porosityCells{numberOfSlices}(timeStep) * ...
    (x(1, :) >= (numberOfSlices - 1) / numberOfSlices & x(1, :) <= 1)';
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
area = surfaceAreaCells{findSlice(x, numberOfSlices, lengthXAxis)}(timeStep);

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
