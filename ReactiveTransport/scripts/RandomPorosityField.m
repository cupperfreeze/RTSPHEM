%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% General setting variables

% Run dissolution simulation on a 6x3 random porosity field and test
% adaptivity scheme
% cf. [4] Section 5.2

GeometrySizes = sqrt((1-[0.4692, 0.5808, 0.7421, 0.7212, 0.4265, 0.7947, 0.6197, 0.7475, 0.376, ...
    0.516, 0.6151, 0.6279, 0.7846, 0.5311, 0.7902, 0.7807, 0.6057, 0.7722]')./4);
global EPS;
EPS = eps;
global Solver
Solver = 'BlockPrec';

numberOfSlicesX = 6;
numberOfSlicesY = 3;
numberOfSlices = numberOfSlicesX * numberOfSlicesY;
multi = 1;
lengthXAxis = 1 * 0.1; % [ cm ]
lengthYAxis = 0.5 * 0.1; % [ cm ]
numPartitionsMicroscale = 64; % Number of partitions in each direction
macroscaleStepSize = lengthXAxis / 100; % [ dm ]
meanGrainSize = 10^(-4); % [ dm ]


tic; % Preprocessing

% Physical parameters
dimension = 2;
spaceScaleFactor = 10 * sqrt(numberOfSlices); % spaceScaleFactor [length of Y] = 1 [cm]

% Parameters in dm and s
dynamicViscosity = 1100e-7; % [ kg dm^(-1) s^(-1) ]

diffusionCoefficient = 0.2 * 1e-4 * 1e-2; % [ dm^2 s^(-1) ]
molarVolume = 64.12 * 1e-3; % [ dm^3 mol^(-1) ]
%rateCoefficientHydrogen = 0.89 * 1e-2; % [ mol dm^(-2) s^(-1) ]
rateCoefficientTST = 4.5e-4 * 1e-2; % [ mol dm^(-2) s^(-1) ]
equilibriumRateConstant = 10^(-16.5);
inletVelocity = 0.01 * 1e-1; % [ dm s^(-1) ]; Pe = 5000
%inletVelocity = 0.001 * 1e-1; % [ dm s^(-1) ]; Pe = 50

initialHydrogenConcentration = 1e-5; % [ mol dm^(-3) ]

dissolutionReactionRate = @(cHydrogen, cCalcium, cCarbonate, cMagnesium) ...
    (cHydrogen^0.5 * rateCoefficientTST) ...
    * (1 - cCalcium * cCarbonate^2 * cMagnesium / equilibriumRateConstant);

% Simulation and discretization parameters (comment for convergence tests)

recalculated = cell(numberOfSlices, 1);
[recalculated{:}] = deal(zeros(numberOfSlices, 1));
disp(['numSlices = ', num2str(numberOfSlices)]);
pecletNumber = (inletVelocity * lengthXAxis) / diffusionCoefficient;
fprintf(['Peclet number: ', num2str(pecletNumber), '\n']);

% dt_macro = 1 * 3600; % Pe = 5000
% initialMacroscaleTimeStepSize = 160; % Pe = 50
initialMacroscaleTimeStepSize = 60 * 2^(-4);
% dt_micro = dt_macro / 100; % 0;
numMicroscaleTimeSteps = 10 / multi;

endTime = 8 * 600 * 3600; % [h] Pe = 5000
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
        %        timeStepSizeFactor = 1.5;
        timeStepSizeFactor = 2^(1 / multi);
        numTimeSteps = floor(log(endTime / initialMacroscaleTimeStepSize) ...
            /log(timeStepSizeFactor));
        timeSteps = [0; timeStepSizeFactor.^(0:numTimeSteps)'] ...
            * initialMacroscaleTimeStepSize;
        if (timeSteps(end) < endTime)
            timeSteps(end+1) = endTime;
        end

    case 'expmax'
        maximalStep = 10^6;
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
%timeSteps = [timeSteps; timeSteps(end)+500000.*[1,2,3,4]']
numTimeSlices = numel(timeSteps);
levelSetEvolutionTime = NaN(numTimeSlices, 1);
cellProblemTime = NaN(numTimeSlices, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization of level set method variables

microscaleGrid = FoldedCartesianGrid(dimension, ...
    kron(ones(1, dimension), [-0.5, 0.5]), ...
    numPartitionsMicroscale*ones(1, dimension));
coord = microscaleGrid.coordinates;
coordCell = mat2cell(coord, ones(1, microscaleGrid.nodes), dimension);
helpGridHyPHM = Grid(coord, microscaleGrid.triangles); %same Grid unfolded in HyPHM format

clear coord;

initialLevelSetDataCells = cell(numberOfSlices, 1);
%lol=[0.1,0.1,0.1,0.1,0.1,0.1, 0.12, 0.12,0.12,0.12,0.12,0.12,0.15, 0.15, 0.15, 0.15, 0.15, 0.15];
for i = 1:numberOfSlices
    temp = GeometrySizes(i) - 0.1; % 0.3*rand();
    initialLevelSetFunc = @(x) 0.1 + temp - norm(x, inf);
    initialLevelSetDataCells{i} = cellfun(initialLevelSetFunc, coordCell);
end
% maybe [ initialLevelSetDataCells{:} ] = deal( cellfun( ...
%           initialLevelSetFunc, coordCell ) ); -> better (however)?

% X = microscaleGrid.reshape( coord(:,1) );
% Y = microscaleGrid.reshape( coord(:,2) );

interfaceNormalVelocity = @(t, x, cHydrogen, cCalcium, cCarbonate, cMagnesium) ...
    spaceScaleFactor * molarVolume * dissolutionReactionRate(cHydrogen, ...
    cCalcium, cCarbonate, cMagnesium);

currentTime = 0;

currentLevelSetDataCells = cell(numberOfSlices, 1);

levelSet = cell(numberOfSlices, 1);
for i = 1:numberOfSlices
    currentLevelSetDataCells{i} = initialLevelSetDataCells{i};
    levelSetData = NaN(numel(initialLevelSetDataCells{i}), 1);
    levelSetData(:, 1) = currentLevelSetDataCells{i};
    levelSet{i} = levelSetData;
end
oldLevelSetDataCells = currentLevelSetDataCells;
clear levelSetData;


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
clear numTimeSlices;

surfaceArea = cell(numberOfSlices, 1);
triangleVolumesOld = cell(numberOfSlices, 1);

for i = 1:numberOfSlices

    [cellProblemSystemMatrix, rhs, isDoF, triangleVolumes, triangleSurfaces] ...
        = assembleCellProblem(microscaleGrid, levelSet{i}(:, 1));

    triangleVolumesOld{i} = triangleVolumes; %reference for adaptive treatment of cell problems
    SOL = solveSystemFE(microscaleGrid, cellProblemSystemMatrix, rhs, isDoF);
    [diffusion, porosities{i}(1)] ...
        = computeDiffusionTensor(microscaleGrid, SOL, triangleVolumes);
    diffusionTensors{i}(:, 1) = diffusion(:);

    % porosities{i}(1) = porosity;
    permeability = computePermeabilityTensor(helpGridHyPHM, levelSet{i}(:, 1));
    permeabilityTensors{i}(:, 1) = permeability(:) * (meanGrainSize^2);
    surfaceArea{i}(1) = sum(triangleSurfaces);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization of HyPHM variables (cont.)

flowStepper = Stepper(0:1);
transportStepper = Stepper(timeSteps);

porosityFunc = @(t, x) porosityHelperFun(t, x, porosities, numberOfSlicesX, numberOfSlicesY, 1, lengthXAxis, lengthYAxis);
diffusionFunc = @(t, x) diffusionHelperFun(t, x, diffusionTensors, ...
    numberOfSlicesX, numberOfSlicesY, 1, lengthXAxis, lengthYAxis) * diffusionCoefficient;
permeabilityFunc = @(t, x) permeabilityHelperFun(t, x, permeabilityTensors, ...
    numberOfSlicesX, numberOfSlicesY, 1, lengthXAxis, lengthYAxis);

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
darcy.computeLevel('s');

temp1 = darcy.Q.getdata(1);
temp2 = temp1;
for ed = 1:gridHyPHM.numE
    porosityFunc(0, mean(gridHyPHM.coordV(gridHyPHM.V0E(ed, :), :)));
    temp2(ed) = temp1(ed); % ./porosityFunc(0,mean(gridHyPHM.coordV(gridHyPHM.V0E(ed,:),:)));
end
flow = Variable(gridHyPHM, flowStepper, 'Flow', 'RT0');
flow.setdata(0, temp2); %variable containing actual flow velocity for normal velocity modification
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
calciumTransport.F = Variable(gridHyPHM, transportStepper, 'F', 'P0');
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
carbonateTransport.F = Variable(gridHyPHM, transportStepper, 'F', 'P0');
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
magnesiumTransport.F = Variable(gridHyPHM, transportStepper, 'F', 'P0');
magnesiumTransport.isUpwind = 'exp';

disp(' ');
disp('Initializing hydrogen transport...');

hydrogenConcentration = Variable(gridHyPHM, transportStepper, ...
    'H^+_pH', 'P0');
hydrogenConcentration.setdata(0, @(t, x) initialHydrogenConcentration*max(inletVelocity - 0.001 * x(1), 0));
hydrogenTransport = Transport(gridHyPHM, transportStepper, 'H^+ Transport');
hydrogenTransport.id2N = {2};
hydrogenTransport.id2F = {1, 3, 4};
%hydrogenTransport.id2D = {2,3,1,4};
hydrogenTransport.U = hydrogenConcentration;
hydrogenTransport.gF.setdata( ...
    @(t, x) -initialHydrogenConcentration*inletVelocity*(x(1) < EPS));
% hydrogenTransport.gF.setdata( @(t,x) -1e-7 );
%hydrogenTransport.uD.setdata( @(t,x) initialHydrogenConcentration );
hydrogenTransport.A.setdata(0, porosityFunc);
hydrogenTransport.C.setdata(0, darcy.Q.getdata(1));
hydrogenTransport.D = Variable(gridHyPHM, transportStepper, 'Diffusion', 'P0P0P0P0');
hydrogenTransport.D.setdata(0, diffusionFunc);
hydrogenTransport.isUpwind = 'exp';
speciesCells = {calciumTransport; carbonateTransport; hydrogenTransport; magnesiumTransport};

calciumDataOld = calciumConcentration.getdata(0);
carbonateDataOld = carbonateConcentration.getdata(0);
magneisumDataOld = magnesiumConcentration.getdata(0);
hydrogenDataOld = hydrogenConcentration.getdata(0);

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

    % average concentrations and Darcy velocity over subdomains related to
    % unit cells
    for i = 1:numberOfSlicesX
        for j = 1:numberOfSlicesY

            isInArea = ...
                (gridHyPHM.baryT(:, 1) >= ((i - 1) / numberOfSlicesX) * lengthXAxis) & ...
                (gridHyPHM.baryT(:, 1) < (i / numberOfSlicesX) * lengthXAxis) & ...
                (gridHyPHM.baryT(:, 2) >= ((j - 1) / numberOfSlicesY) * lengthYAxis) & ...
                (gridHyPHM.baryT(:, 2) < (j / numberOfSlicesY) * lengthXAxis);
            index = (j - 1) * numberOfSlicesX + i;
            calciumSliceAverageCells{index} = sum(calciumData(isInArea).* ...
                gridHyPHM.areaT(isInArea)) / sum(gridHyPHM.areaT(isInArea));
            carbonateSliceAverageCells{index} = sum(carbonateData(isInArea).* ...
                gridHyPHM.areaT(isInArea)) / sum(gridHyPHM.areaT(isInArea));
            magnesiumSliceAverageCells{index} = sum(magnesiumData(isInArea).* ...
                gridHyPHM.areaT(isInArea)) / sum(gridHyPHM.areaT(isInArea));
            hydrogenSliceAverageCells{index} = sum(hydrogenData(isInArea).* ...
                gridHyPHM.areaT(isInArea)) / sum(gridHyPHM.areaT(isInArea));
            flowSliceAverageCells{index} = sum(flowVelocitiesCartesian(isInArea, :).* ...
                gridHyPHM.areaT(isInArea), 1) / sum(gridHyPHM.areaT(isInArea));
        end
    end


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

        microscaleTimeStepSize = transportStepper.curtau / numMicroscaleTimeSteps;
        microscaleTime = currentTime;
        oldMicroscaleTime = currentTime;
        for j = 1:numMicroscaleTimeSteps
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

    %     t_old = currentTime;
    %     currentTime = currentTime + dt_macro;
    %     oldLevelSetDataCells = currentLevelSetDataCells;
    currentTime = currentTime + macroscaleTimeStepSize;

    levelSetEvolutionTime(timeIterationStep) = toc;
    disp(['    ... done in ', ...
        num2str(levelSetEvolutionTime(timeIterationStep)), ' seconds.']);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Calculation of effective parameters
    %if (mod(timeIterationStep,multi) == 0 || abs(timeIterationStep - length(timeSteps) + 1)<eps )    %every multi-th and last iteration
    if true
        disp('Calculation of effective parameters ...');
        tic;
        count = 0;
        loccount = zeros(numberOfSlices, 1);

        for i = 1:numberOfSlices

            [cellProblemSystemMatrix, rhs, isDoF, triangleVolumes, ...
                triangleSurfaces] = assembleCellProblem(microscaleGrid, ...
                levelSet{i}(:, timeIterationStep + 1));

            if sum(abs(triangleVolumesOld{i} - triangleVolumes) > 0.1*2/numPartitionsMicroscale^2) > eps % if recalculation necessary
                loccount(i) = 1;

                SOL = solveSystemFE(microscaleGrid, cellProblemSystemMatrix, rhs, isDoF);
                [diffusion, porosities{i}(timeIterationStep + 1)] ...
                    = computeDiffusionTensor(microscaleGrid, SOL, ...
                    triangleVolumes);
                diffusionTensors{i}(:, timeIterationStep + 1) = diffusion(:);
                %             porosities{i}( timeIterationStep ) = porosity;

                permeability = computePermeabilityTensor(helpGridHyPHM, ...
                    levelSet{i}(:, timeIterationStep + 1));
                permeabilityTensors{i}(:, timeIterationStep + 1) = permeability(:) * (meanGrainSize^2);
            else % copy old data
                diffusionTensors{i}(:, timeIterationStep + 1) = diffusionTensors{i}(:, timeIterationStep);
                porosities{i}(timeIterationStep + 1) = porosities{i}(timeIterationStep);
                permeabilityTensors{i}(:, timeIterationStep + 1) = permeabilityTensors{i}(:, timeIterationStep);
                count = count + 1;
            end
            surfaceArea{i}(timeIterationStep + 1) = sum(triangleSurfaces);

        end
        recalculated{timeIterationStep, 1}(:) = loccount;
        disp([num2str(count), ' recalculation(s) saved by adaptivity :-)']);
        cellProblemTime(timeIterationStep) = toc;
        disp(['    ... done in ', num2str(cellProblemTime(timeIterationStep)), ...
            ' seconds.']);

        permeabilityFunc = @(t, x) permeabilityHelperFun(t, x, permeabilityTensors, ...
            numberOfSlicesX, numberOfSlicesY, timeIterationStep+1, lengthXAxis, lengthYAxis);
        porosityFunc = @(t, x) porosityHelperFun(t, x, porosities, numberOfSlicesX, numberOfSlicesY, timeIterationStep+1, ...
            lengthXAxis, lengthYAxis);
        darcy.D.setdata(1,@(t, x) permeabilityFunc(t, x)./dynamicViscosity);
        flowStepper.next;
        darcy.computeLevel('s');

        temp1 = darcy.Q.getdata(1);

        flow.setdata(0, temp1);
        flowStepper = Stepper(0:1);


    else
        for i = 1:numberOfSlices
            surfaceArea{i}(timeIterationStep + 1) = surfaceArea{i}(timeIterationStep);
            diffusionTensors{i}(:, timeIterationStep + 1) = diffusionTensors{i}(:, timeIterationStep);
            porosities{i}(timeIterationStep + 1) = porosities{i}(timeIterationStep);
            permeabilityTensors{i}(:, timeIterationStep + 1) = permeabilityTensors{i}(:, timeIterationStep);
        end
    end

    magnesiumTransport.C.setdata(timeIterationStep, darcy.Q.getdata(1));
    hydrogenTransport.C.setdata(timeIterationStep, darcy.Q.getdata(1));
    calciumTransport.C.setdata(timeIterationStep, darcy.Q.getdata(1));
    carbonateTransport.C.setdata(timeIterationStep, darcy.Q.getdata(1));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Macroscopic transport step

    porosityFunc = @(t, x) porosityHelperFun(t, x, porosities, numberOfSlicesX, numberOfSlicesY, timeIterationStep+1, lengthXAxis, lengthYAxis);
    diffusionFunc = @(t, x) diffusionHelperFun(t, x, diffusionTensors, ...
        numberOfSlicesX, numberOfSlicesY, timeIterationStep+1, lengthXAxis, lengthYAxis) * diffusionCoefficient;
    permeabilityFunc = @(t, x) permeabilityHelperFun(t, x, permeabilityTensors, ...
        numberOfSlicesX, numberOfSlicesY, timeIterationStep+1, lengthXAxis, lengthYAxis);

    %     surfaceAreaFun = @(x) areaHelperFun( -1, x, microscaleGrid, levelSet, numberOfSlices, timeStep );
    surfaceAreaFunc = @(x) areaHelperFun(-1, x, surfaceArea, numberOfSlicesX, numberOfSlicesY, ...
        timeIterationStep+1, lengthXAxis, lengthYAxis);

    specificSurfaceArea = cellfun(surfaceAreaFunc, macroCoordCell) * spaceScaleFactor;
    % TODO evaluate dissolutionReactionRate and other related function
    % handles (really needed?)
    calciumTransportRhsData = 0 * specificSurfaceArea ...
        .* (sqrt(hydrogenData) * rateCoefficientTST) .* ...
        (1 - calciumData .* carbonateData .* carbonateData .* magnesiumData ./ equilibriumRateConstant);

    carbonateTransportRhsData = 2 * calciumTransportRhsData;
    magnesiumTransportRhsData = calciumTransportRhsData;

    hydrogenTransportRhsData = 0 * calciumTransportRhsData;

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
    %     surfaceArea = surfaceIntegral( g, levelSet{1}( :, timeStep ), @(x) 1 ) / ...
    %         spaceScaleFactor;

    %     calciumRHS = 47.50 * k_3 * ( 1 - calciumData .* carbonateData / K_s );
    %     calciumRHS = surfaceArea * k_3 * ( 1 - calciumData .* carbonateData / K_s );
    %     carbonateRHS = calciumRHS;

    calciumTransport.A.setdata(timeIterationStep, porosityFunc);
    calciumTransport.D.setdata(timeIterationStep, diffusionFunc);
    calciumTransport.F.setdata(timeIterationStep, calciumTransportRhsData);
    carbonateTransport.A.setdata(timeIterationStep, porosityFunc);
    carbonateTransport.D.setdata(timeIterationStep, diffusionFunc);
    carbonateTransport.F.setdata(timeIterationStep, carbonateTransportRhsData);
    hydrogenTransport.A.setdata(timeIterationStep, porosityFunc);
    hydrogenTransport.D.setdata(timeIterationStep, diffusionFunc);
    hydrogenTransport.F.setdata(timeIterationStep, hydrogenTransportRhsData);
    magnesiumTransport.A.setdata(timeIterationStep, porosityFunc);
    magnesiumTransport.D.setdata(timeIterationStep, diffusionFunc);
    magnesiumTransport.F.setdata(timeIterationStep, magnesiumTransportRhsData);

    newtonIteration(speciesCells, nonlinearFunc, nonlinearJacFunc, 10);

end % while

calciumConcentration.visualize;
carbonateConcentration.visualize;
hydrogenConcentration.visualize;
darcy.U.visualize();

% dataFolder = 'data/Molins2017/ConvergenceTests_201805/MacroGridConvergence/Newton/';

% filePrefix = [ dataFolder, 'Physics_Data_', ...
%     '_partsMacro_Peclet_', num2str( pecletNumber ) ];

dateString = datestr(now, 'yyyymmdd_HHMM');
% try
%     save( [ dataFolder, 'Physics_Data_', dateString, '.mat' ], ...
%         'dimension', 'spaceScaleFactor', 'diffusionCoefficient', 'molarVolume', ...
%         'rateCoefficientHydrogen', 'rateCoefficientTST', 'equilibriumRateConstant', ...
%         'inletVelocity', 'initialHydrogenConcentration', 'dissolutionReactionRate', ...
%         'lengthXAxis', 'lengthYAxis', 'pecletNumber', 'endTime', '-v7.3' );
% catch
%     disp( 'Physics data could not be saved.' );
% end
% try
%     save( [ dataFolder, 'Level_set_', dateString, '.mat' ], ...
%         'initialLevelSetFunc', 'interfaceNormalVelocity', 'levelSet', '-v7.3' );
% catch
%     disp( 'Level set data could not be saved.' );
% end
% try
%     save( [ dataFolder, 'Simulation_data_', dateString, '.mat' ], ...
%         'numPartitionsMicroscale', 'macroscaleStepSize', 'numberOfSlices', ...
%         'pecletNumber', 'initialMacroscaleTimeStepSize', 'numMicroscaleTimeSteps', ...
%         'timeStepperType', 'timeSteps', 'timeStepSizeFactor', 'microscaleGrid', ...
%         'gridHyPHM', 'flowStepper', 'transportStepper', '-v7.3' );
% catch
%     disp( 'Simulation data could not be saved.' );
% end
% try
%     save( [ dataFolder, 'General_results_', dateString, '.mat' ], ...
%         'diffusionTensors', 'porosities', 'surfaceArea', 'flow', ...
%         'levelSetEvolutionTime', 'cellProblemTime', 'preprocessingTime', '-v7.3' );
% catch
%     disp( 'General results could not be saved.' );
% end
% try
%     save( [ dataFolder, 'Calcium_', dateString, '.mat' ], ...
%         'calciumConcentration', 'calciumTransport', '-v7.3' );
% catch
%     disp( 'Calcium data could not be saved.' );
% end
% try
%     save( [ dataFolder, 'Carbonate_', dateString, '.mat' ], ...
%         'carbonateConcentration', 'carbonateTransport', '-v7.3' );
% catch
%     disp( 'Carbonate data could not be saved.' );
% end
% try
%     save( [ dataFolder, 'Hydrogen_', dateString, '.mat' ], ...
%         'hydrogenConcentration', 'hydrogenTransport', '-v7.3' );
% catch
%     disp( 'Hydrogen data could not be saved.' );
% end
%
% diary off;
% movefile( logFilePath, strcat( dataFolder, 'Log_', dateString ) );

% save( [ filePrefix, '_Calcium_', dateString, '.mat' ], ...
%     'calciumConcentration', 'calciumData', 'calciumDataOld', ...
%     'calciumSliceAverageCells', 'calciumTransport', ...
%     'calciumTransportRhsData' );
% save( [ filePrefix, '_Carbonate_', dateString, '.mat' ], ...
%     'carbonateConcentration', 'carbonateData', 'carbonateDataOld', ...
%     'carbonateSliceAverageCells', 'carbonateTransport', ...
%     'carbonateTransportRhsData' );
% save( [ filePrefix, '_Hydrogen_', dateString, '.mat' ], ...
%     'hydrogenConcentration', 'hydrogenData', 'hydrogenDataOld', ...
%     'hydrogenSliceAverageCells', 'hydrogenTransport', ...
%     'hydrogenTransportRhsData' );
% save( [ dataFolder, 'Data_pH_', num2str( 1 / macroscaleStepSize ), ...
%     '_partsMacro_LevelSet_Peclet_50_', dateString, '.mat' ], ...
%     'levelSet', '-v7.3' );
% save( [ filePrefix, '_Other_', dateString, '.mat' ], ...
%     'cellProblemSystemMatrix', 'cellProblemTime', 'coord', 'coordCell', ...
%     'diffusion', 'diffusionCoefficient', 'diffusionFun', ...
%     'diffusionTensors', 'diffusionTensorsData', 'dimension', ...
%     'dissolutionReactionRate', 'dt_macro', 'dt_micro', 'endTime', ...
%     'equilibriumRateConstant', 'flow', 'flowStepper', 'gridHyPHM', ...
%     'initialData', 'initialHydrogenConcentration', 'inletVelocity', ...
%     'interfaceNormalVelocity', 'isDoF', 'isInArea', 'lengthXAxis', ...
%     'lengthYAxis', 'levelSetData', 'levelSetEvolutionTime', ...
%     'lsf0', 'macroCoordCell', 'macroscaleStepSize', 'microscaleGrid', ...
%     'molarVolume', 'normalSpeed', 'numberOfSlices', ...
%     'numPartitionsMicroscale', 'numTimeSlices', 'phi', 'phiOld', ...
%     'porosities', 'porositiesData', 'porosity', 'porosityFun', ...
%     'preprocessingTime', 'rateCoefficientHydrogen', ...
%     'rateCoefficientTST', 'rateFunction', 'rhs', 'SOL', ...
%     'spaceScaleFactor', 'specificSurfaceArea', 'step', 'surfaceArea', ...
%     'surfaceAreaFun', 't_', 't_micro', 't_old', 'timeCell', 'timeStep', ...
%     'timeSteps', 'transportStepper', 'triangleSurfaces', ...
%     'triangleVolumes', 'velocity', 'X', 'Y' );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Helper functions

function por = porosityHelperFun(t, x, porosityCells, numberOfSlicesX, numberOfSlicesY, timeStep, lengthXAxis, lengthYAxis)
%return porosity for each macroscopic triangle
xSlice = findSlice(x, numberOfSlicesX, lengthXAxis);
ySlice = findSlice([x(2), x(1)], numberOfSlicesY, lengthYAxis);
por = porosityCells{(ySlice-1)*numberOfSlicesX+xSlice}(timeStep);

end


function perm = permeabilityHelperFun(t, x, permeabilityCells, numberOfSlicesX, numberOfSlicesY, timeStep, lengthXAxis, lengthYAxis)
%return permeability tensor for each macroscopic triangle
xSlice = findSlice(x, numberOfSlicesX, lengthXAxis);
ySlice = findSlice([x(2), x(1)], numberOfSlicesY, lengthYAxis);
perm = reshape(permeabilityCells{(ySlice-1) * numberOfSlicesX + xSlice}(:, timeStep), 2, 2);

end

function diff = diffusionHelperFun(t, x, diffusionCells, numberOfSlicesX, numberOfSlicesY, timeStep, lengthXAxis, lengthYAxis)
%return diffusion tensor for each macroscopic triangle
xSlice = findSlice(x, numberOfSlicesX, lengthXAxis);
ySlice = findSlice([x(2), x(1)], numberOfSlicesY, lengthYAxis);
diff = reshape(diffusionCells{(ySlice-1) * numberOfSlicesX + xSlice}(:, timeStep), 2, 2);

end

function area = areaHelperFun(t, x, surfaceAreaCells, numberOfSlicesX, numberOfSlicesY, timeStep, lengthXAxis, lengthYAxis)
%return surface area for each macroscopic triangle
xSlice = findSlice(x, numberOfSlicesX, lengthXAxis);
ySlice = findSlice([x(2), x(1)], numberOfSlicesY, lengthYAxis);
area = surfaceAreaCells{(ySlice-1)*numberOfSlicesX+xSlice}(timeStep);

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


function sliceNumber = findSlice(x, numberOfSlices, lengthAxis)
%return number of slice related to macroscopic coordinate x
sliceNumber = 1 + fix(x(1)./lengthAxis*numberOfSlices);
if sliceNumber > numberOfSlices
    sliceNumber = numberOfSlices;
end
end


function y = normalVelocityModifierFunc(x)

positiveLimitValue = 0.75;
%expFactor = 5e+3; % Value for convergence tests
expFactor = 5e+2;
positive = (x > 0.0);
y = exp(-expFactor*x.^2);
y(positive) = y(positive) + positiveLimitValue ...
    * (1.0 - exp(-expFactor * x(positive).^2));

end
