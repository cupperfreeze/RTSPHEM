% Test of TransportLEVEL and StokesLEVEL against first problem in the benchmark of
%'Simulation of mineral dissolution at the pore scale with evolving
%fluid-solid interfaces: review of approaches and benchmark problem set',
%Molins et al.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% General setting variables
global Solver
Solver = 'StandardDirect';
global EPS;
EPS = eps;

% all quantities in cm
lengthXAxis = 0.1; % [ cm ]
lengthYAxis = 0.05; % [ cm ]
numPartitionsMicroscale = 128; % Number of partitions
macroscaleStepSize = lengthXAxis / 5; % [ dm ]


tic; % Preprocessing

% Physical parameters
dimension = 2;
spaceScaleFactor = 1; % spaceScaleFactor [length of Y] = 1 [cm]

% Parameters in cm and s
diffusionCoefficient = 1e-5; % [ cm^2 s^(-1) ]


rateCoefficientTST = 10^(-4.05); % [ mol cm^(-2) s^(-1) ]

inletVelocity = 0.12; % [ cm s^(-1) ];
initialHydrogenConcentration = 1e-5; % [ mol cm^(-3) ]

% Simulation and discretization parameters (comment for convergence tests)

numberOfSlices = 1;
disp(['numSlices = ', num2str(numberOfSlices)]);
pecletNumber = (inletVelocity * lengthXAxis) / diffusionCoefficient;
fprintf(['Peclet number: ', num2str(pecletNumber), '\n']);

% dt_macro = 1 * 3600; % Pe = 5000
% initialMacroscaleTimeStepSize = 160; % Pe = 50
initialMacroscaleTimeStepSize = 0.002;
% dt_micro = dt_macro / 100; % 0;
numMicroscaleTimeSteps = 10;
endTime = 2; % [h] Pe = 5000
% endTime = 10 * 60 * 60; % Pe = 50
% endTime = 200 * dt_macro;

%% Computation of time steps in simulation
timeStepperType = 'exp';
%timeStepperType = 'expmax';

switch (timeStepperType)
    case 'linear'
        if (mod(endTime, initialMacroscaleTimeStepSize) < EPS)
            timeSteps = 0:initialMacroscaleTimeStepSize:endTime;
        else
            timeSteps = [0:initialMacroscaleTimeStepSize:endTime, endTime];
        end
        timeStepSizeFactor = 1;
    case 'exp'
        timeStepSizeFactor = 2^(1 / 3);
        %  timeStepSizeFactor = 2;
        numTimeSteps = floor(log(endTime / initialMacroscaleTimeStepSize) ...
            /log(timeStepSizeFactor));
        timeSteps = ([0; (timeStepSizeFactor.^(0:numTimeSteps)')] ...
            * initialMacroscaleTimeStepSize);
        if (timeSteps(end) < endTime)
            timeSteps(end+1) = endTime;
        end
        size(timeSteps)
    case 'expmax'
        maximalStep = 5 * 10^5;
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

%timeSteps = [(0:initialMacroscaleTimeStepSize:0.2),((0.2 + initialMacroscaleTimeStepSize):5*initialMacroscaleTimeStepSize:3) ]
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


initialLevelSetFunc = @(x) 0.01 - norm(x-[0.05, 0.025]);
initialLevelSetDataCells = cell(numberOfSlices, 1);
initialLevelSetDataCells{1} = cellfun(initialLevelSetFunc, coordCell);

[a, b] = meshgrid(0:lengthYAxis/numPartitionsMicroscale:lengthXAxis, 0:lengthYAxis/numPartitionsMicroscale:lengthYAxis);
contour(a, b, reshape(initialLevelSetDataCells{1}, 2 * numPartitionsMicroscale + 1, numPartitionsMicroscale + 1)', [0, 1])
axis equal


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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gridHyPHM = Grid(coord, microscaleGrid.triangles);
% [ gridHyPHM,   levelSet{1}] = localMeshRefinement(gridHyPHM, ...
%     levelSet{1}(:,1) );

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

for i = 1:numberOfSlices

    [cellProblemSystemMatrix, rhs, isDoF, triangleVolumes, triangleSurfaces] ...
        = assembleCellProblem(microscaleGrid, levelSet{i}(:, 1));

    SOL = solveSystemFE(microscaleGrid, cellProblemSystemMatrix, rhs, isDoF);
    [diffusion, porosities{i}(1)] ...
        = computeDiffusionTensor(microscaleGrid, SOL, ...
        triangleVolumes);
    surfaceArea{i}(1) = sum(triangleSurfaces);
    Volume = (0.05 * 0.1 - porosities{1}(1))

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


flow = Variable(gridHyPHM, flowStepper, 'Flow', 'RT0');
helper = StokesL.U.getdata(1);
flow.setdata(helper(gridHyPHM.numV + 1:end, 1).*gridHyPHM.nuE(:, 1)+helper(gridHyPHM.numV + 1:end, 2).*gridHyPHM.nuE(:, 2));
porosityFunc = @(t, x) porosityHelperFun(t, x, porosities, numberOfSlices, 1);


% Boundary IDs: 1 = down, 2 = right, 3 = up, 4 = left

disp(' ');
disp('Initializing calcium transport...');


calciumConcentration = Variable(gridHyPHM, transportStepper, ...
    'Ca^(2+)', 'P0');
calciumConcentration.setdata(0, @(t, x) 0);
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
calciumTransport.gF.setdata(0, zeros(gridHyPHM.numE, 1));
calciumTransport.A = Variable(gridHyPHM, transportStepper, 'Porosity', 'P0');
calciumTransport.A.setdata(0, @(t, x) 1);
% calciumTransport.C.setdata( P2P2.P2P2toRT0slice( g_HyPHM, flow.getdata(1) ) );
calciumTransport.C.setdata(0, flow.getdata(1));
calciumTransport.isUpwind = 'exp';

hydrogenConcentration = Variable(gridHyPHM, transportStepper, ...
    'H^+_pH', 'P0');
hydrogenConcentration.setdata(0, @(t, x) initialHydrogenConcentration);
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
hydrogenTransport.A.setdata(0, @(t, x) 1);
hydrogenTransport.C.setdata(0, flow.getdata(1));
hydrogenTransport.isUpwind = 'exp';


calciumDataOld = calciumConcentration.getdata(0);
hydrogenDataOld = hydrogenConcentration.getdata(0);

preprocessingTime = toc; % Preprocessing

disp(' ');
disp(['Preprocessing done in ', num2str(preprocessingTime), ' seconds.']);

hydrogenTransport.L.setdata(0, [levelSet{1}(:, 1)]);
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
    if min(hydrogenTransport.U.getdata(timeIterationStep - 1)) < -eps
        disp('negative H+ concentration detected')
        break
    end
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

    currentTime = currentTime + macroscaleTimeStepSize;


    % Level set evolution step


    levelSet{1}(:, timeIterationStep + 1) = levelSet{1}(:, timeIterationStep);

    hydrogenTransport.L.setdata(timeIterationStep, levelSet{1}(:, timeIterationStep + 1));
    calciumTransport.L.setdata(timeIterationStep, levelSet{1}(:, timeIterationStep + 1));


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Calculation of effective parameters

    disp('Calculation of effective parameters ...');
    tic;
    %
    %     for i = 1:numberOfSlices
    %
    %         [ cellProblemSystemMatrix, rhs, isDoF, triangleVolumes, ...
    %             triangleSurfaces ] = assembleCellProblem( microscaleGrid, ...
    %             levelSet{i}( :, timeIterationStep ) );
    %
    %         surfaceArea{i}( timeIterationStep ) = sum( triangleSurfaces );
    %
    %     end
    cellProblemTime(timeIterationStep) = toc;
    disp(['    ... done in ', num2str(cellProblemTime(timeIterationStep)), ...
        ' seconds.']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Stokes
    %Calculate Pressure distribution and Stokes velocity
    disp(' ');
    disp('Initializing Stokes problem...');


    StokesL.L.setdata(levelSet{1}(:, timeIterationStep + 1))

    flowStepper.next;

    helper = StokesL.U.getdata(1);
    flow.setdata(helper(gridHyPHM.numV + 1:end, 1).*gridHyPHM.nuE(:, 1)+helper(gridHyPHM.numV + 1:end, 2).*gridHyPHM.nuE(:, 2));
    flowStepper.prev;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Macroscopic transport step


    %     surfaceAreaFun = @(x) areaHelperFun( -1, x, microscaleGrid, levelSet, numberOfSlices, timeStep );
    surfaceAreaFunc = @(x) areaHelperFun(-1, x, surfaceArea, numberOfSlices, ...
        timeIterationStep);

    %    specificSurfaceArea = cellfun( surfaceAreaFunc, macroCoordCell ) * spaceScaleFactor;
    % TODO evaluate dissolutionReactionRate and other related function
    % handles (really needed?)
    Levels = hydrogenTransport.L.getdata(timeIterationStep);
    %%%%%%
    % handle edge orientation
    edgeOrientation = zeros(gridHyPHM.numE, 1);
    for kt = 1:gridHyPHM.numT
        if sum(Levels(gridHyPHM.V0T(kt, :)) > -eps) > 2.5 %solid triangle
            edgeOrientation(gridHyPHM.E0T(kt, :)) = gridHyPHM.sigE0T(kt, :);
        end
    end


    %%%%%%


    %concentrations P_0(T)--> P_0(E)
    calciumDataE = zeros(gridHyPHM.numE, 1);
    hydrogenDataE = zeros(gridHyPHM.numE, 1);


    calciumData = calciumConcentration.getdata(timeIterationStep-1);
    hydrogenData = hydrogenConcentration.getdata(timeIterationStep-1);

    for i = 1:gridHyPHM.numE
        if ((Levels(gridHyPHM.V0E(i, 1)) > -eps) & (Levels(gridHyPHM.V0E(i, 2)) > -eps))
            helper = gridHyPHM.T0E(i, :);
            calciumDataE(i, 1) = max(calciumData(helper));
            hydrogenDataE(i, 1) = max(hydrogenData(helper));
        end
    end


    %%%%%%


    calciumTransport.A.setdata(timeIterationStep, @(t, x) 1);
    calciumTransport.gF.setdata(timeIterationStep, zeros(gridHyPHM.numE, 1));
    calciumTransport.C.setdata(timeIterationStep, flow.getdata(1));

    hydrogenTransport.A.setdata(timeIterationStep, @(t, x) 1);
    hydrogenTransport.gF.setdata(timeIterationStep, hydrogenTransport.gF.getdata(0));
    hydrogenTransport.C.setdata(timeIterationStep, flow.getdata(1));

    SorceScale = zeros(gridHyPHM.numT, 1);
    levels = levelSet{1}(:, 1);
    for kT = 1:gridHyPHM.numT
        L = 0;
        for i = 1:3
            if levels(gridHyPHM.V0E(gridHyPHM.E0T(kT, i), 1)) > -eps & levels(gridHyPHM.V0E(gridHyPHM.E0T(kT, i), 2)) > -eps
                L = L + gridHyPHM.areaE(gridHyPHM.E0T(kT, i)) / gridHyPHM.areaT(kT);
            end
        end
        SorceScale(kT) = L;
    end

    speciesCells = {hydrogenTransport, calciumTransport};
    nonlinearFunc = cell(2, 1);
    nonlinearFunc{1} = @(x, y) -SorceScale ...
        .* (x * rateCoefficientTST * 1000);
    nonlinearFunc{2} = @(x, y) -nonlinearFunc{1}(x, y);

    nonlinearJacFunc = cell(2, 2);
    nonlinearJacFunc{1, 1} = @(x, y) -SorceScale ...
        .* (rateCoefficientTST * 1000);
    nonlinearJacFunc{1, 2} = @(x, y) SorceScale * 0;

    nonlinearJacFunc{2, 1} = @(x, y) -nonlinearJacFunc{1, 1}(x, y);
    nonlinearJacFunc{2, 2} = @(x, y) nonlinearJacFunc{1, 2}(x, y);
    newtonIteration(speciesCells, nonlinearFunc, nonlinearJacFunc, 2);


    Rate(timeIterationStep) = (calciumTransport.Q.getdata(timeIterationStep) .* (gridHyPHM.baryE(:, 1) > (0.1 - eps)))' * gridHyPHM.areaE / surfaceArea{1};
    hydrogenData = hydrogenConcentration.getdata(timeIterationStep);
    for i = 1:gridHyPHM.numE
        if gridHyPHM.baryE(i, 1) > 0.1 - eps
            helper = max(gridHyPHM.T0E(i, :));
            hydrogenDataE(i, 1) = hydrogenData(helper);
        end
    end
    cOut(timeIterationStep) = (hydrogenDataE .* (gridHyPHM.baryE(:, 1) > (0.1 - eps)))' * gridHyPHM.areaE / 0.05;

    % hydrogenTransport.U.setdata(timeIterationStep, max(hydrogenTransport.U.getdata(timeIterationStep),0));

end % while

% plot results related to figure 5
figure
plot(timeSteps(2:end), Rate);
xlabel('time [s]')
ylabel('average dissolution rate')

figure
plot(timeSteps(2:end), cOut);
xlabel('time [s]')
ylabel('effluent concentration')

%% Helper functions

function por = porosityHelperFun(t, x, porosityCells, numberOfSlices, timeStep)

por = 0;
for i = 1:numberOfSlices
    por = por + porosityCells{1}(timeStep);
end

%     por = porosityCells{ findSlice( x, numberOfSlices ) }( timeStep );

end

function diff = diffusionHelperFun(t, x, diffusionCells, numberOfSlices, timeStep)

diff = zeros(2);
for i = 1:numberOfSlices
    diff = diff + reshape(diffusionCells{1}(:, timeStep), 2, 2);
end

%     diff = reshape( diffusionCells{ findSlice( x, numberOfSlices ) }( :, ...
%         timeStep ), 2, 2 );

end

function area = areaHelperFun(t, x, surfaceAreaCells, numberOfSlices, timeStep)

area = surfaceAreaCells{findSlice(x, numberOfSlices)}(timeStep);

end


function sliceNumber = findSlice(x, numberOfSlices)

sliceNumber = 1;

end

function timeStep = findTimeStep(t, deltaTime)

timeStep = 1 + fix(t/deltaTime);

end
