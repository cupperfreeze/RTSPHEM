% Test of transport/geometry coupling against second problem in the benchmark of
%'Simulation of mineral dissolution at the pore scale with evolving
%fluid-solid interfaces: review of approaches and benchmark problem set',
%Molins et al.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% General setting variables
global Solver
Solver = 'StandardDirect';
global EPS;
EPS = eps;


lengthXAxis = 0.1; % [ cm ]
lengthYAxis = 0.05; % [ cm ]
numPartitionsMicroscale = 2 * 64; % Number of partitions

tic; % Preprocessing

% Physical parameters
dimension = 2;
spaceScaleFactor = 1; % spaceScaleFactor [length of Y] = 1 [cm]

% Parameters in dm and s
diffusionCoefficient = 1e-5; % [ cm^2 s^(-1) ]
molarVolume = 36.9; % [ cm^3 mol^(-1) ]
%rateCoefficientHydrogen = 0.89 * 1e-2; % [ mol dm^(-2) s^(-1) ]
rateCoefficientTST = 10^(-4.05); % [ mol dm^(-2) s^(-1) ]

%inletVelocity = 0.1* 1e-1; % [ dm s^(-1) ]; Pe = 5000
%inletVelocity = 0.001 * 1e-1; % [ dm s^(-1) ]; Pe = 50
inletVelocity = 0.12; % [ cm s^(-1) ]; Pe = 500
initialHydrogenConcentration = 1e-5; % [ mol cm^(-3) ]

dissolutionReactionRate = @(cHydrogen) ...
    (cHydrogen * 1000 * rateCoefficientTST);

% Simulation and discretization parameters (comment for convergence tests)

numberOfSlices = 1;
disp(['numSlices = ', num2str(numberOfSlices)]);
pecletNumber = (inletVelocity * lengthXAxis) / diffusionCoefficient;
fprintf(['Peclet number: ', num2str(pecletNumber), '\n']);

% dt_macro = 1 * 3600; % Pe = 5000
% initialMacroscaleTimeStepSize = 160; % Pe = 50
initialMacroscaleTimeStepSize = 1;
% dt_micro = dt_macro / 100; % 0;
numMicroscaleTimeSteps = 1;
endTime = 2700; % [h] Pe = 5000
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
        %  timeStepSizeFactor = sqrt(2);
        timeStepSizeFactor = 2;
        numTimeSteps = floor(log(endTime / initialMacroscaleTimeStepSize) ...
            /log(timeStepSizeFactor));
        timeSteps = [0; timeStepSizeFactor.^(0:numTimeSteps)'] ...
            * initialMacroscaleTimeStepSize;
        if (timeSteps(end) < endTime)
            timeSteps(end+1) = endTime;
        end
    case 'expmax'
        maximalStep = 90;
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


initialLevelSetFunc = @(x) 0.01 - norm(x-[0.05, 0.025]);
initialLevelSetDataCells = cell(numberOfSlices, 1);
initialLevelSetDataCells{1} = cellfun(initialLevelSetFunc, coordCell);


[a, b] = meshgrid(0:lengthYAxis/numPartitionsMicroscale:lengthXAxis, 0:lengthYAxis/numPartitionsMicroscale:lengthYAxis);
contour(a, b, reshape(initialLevelSetDataCells{1}, 2 * numPartitionsMicroscale + 1, numPartitionsMicroscale + 1)', [0, 1])
axis equal

interfaceNormalVelocity = @(cHydrogen) ...
    molarVolume * dissolutionReactionRate(cHydrogen);

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


% assemble HyPHM grid and label outer edges
gridHyPHM = Grid(coord, microscaleGrid.triangles);


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
%PorositiyEstimate = sum(levelSet{1}(:,1)<-eps)/numel(levelSet{1}(:,1))

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


flow = Variable(gridHyPHM, flowStepper, 'Flow', 'RT0');
helper = StokesL.U.getdata(1);
flow.setdata(helper(gridHyPHM.numV + 1:end, 1).*gridHyPHM.nuE(:, 1)+helper(gridHyPHM.numV + 1:end, 2).*gridHyPHM.nuE(:, 2));
porosityFunc = @(t, x) porosityHelperFun(t, x, porosities, numberOfSlices, 1);


% Boundary IDs: 1 = down, 2 = right, 3 = up, 4 = left


disp(' ');
disp('Initializing hydrogen transport...');

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

hydrogenDataOld = hydrogenConcentration.getdata(0);

preprocessingTime = toc; % Preprocessing

disp(' ');
disp(['Preprocessing done in ', num2str(preprocessingTime), ' seconds.']);

hydrogenTransport.L.setdata(0, [levelSet{1}(:, 1)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mapping Vertices --> traingles

VertexTriMatrix = cell(1, gridHyPHM.numV);
for i = 1:gridHyPHM.numT
    for j = 1:3
        VertexTriMatrix{1, gridHyPHM.V0T(i, j)} = [VertexTriMatrix{1, gridHyPHM.V0T(i, j)}, i];
    end
end
Rate = nan(1, numel(timeSteps));
Volume = nan(1, numel(timeSteps));

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


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   inner corrector
    if timeIterationStep > 0
        speciesCells = {hydrogenTransport};
        [out, Corrector] = Continuation(speciesCells, timeIterationStep, levelSet{1}(:, timeIterationStep), 20);
        hydrogenTransport = out{1};
    else
        Corrector = zeros(gridHyPHM.numT, 1);
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Level set evolution step

    disp('Evolution of level set ...');
    tic;
    coordCell = mat2cell(microscaleGrid.coordinates, ...
        ones(1, microscaleGrid.nodes), dimension);
    timeCell = num2cell(currentTime*ones(microscaleGrid.nodes, 1));


    hydrogenData = hydrogenTransport.U.getdata(timeIterationStep-1);


    %Convert Concentrations from T --> V

    hydrogenDataV = zeros(gridHyPHM.numV, 1);

    for i = 1:gridHyPHM.numV
        %[row,~] = find(gridHyPHM.V0T(:,:)==i);
        row = VertexTriMatrix{1, i};

        hydrogenDataV(i) = min(hydrogenData(row));
    end
    %
    levels = levelSet{1}(:, timeIterationStep);
    %     for i=1:gridHyPHM.numT
    %        if sum(levels(gridHyPHM.V0T(i,:))>-eps)>1.5 & sum(levels(gridHyPHM.V0T(i,:))>-eps)<2.5
    %             hydrogenDataV(gridHyPHM.V0T(i,:)) = hydrogenData(i)  ;
    %        end
    %     end

    normalSpeed = arrayfun(interfaceNormalVelocity, ...
        hydrogenDataV);

    normalSpeedMax = max(abs(normalSpeed));

    %     InterfacePoints = getInterfacePoints(levelSet{1}(:,timeIterationStep));
    %     normalSpeednew = zeros(numel(levelSet{1}(:,timeIterationStep)),numel(InterfacePoints));
    %     for i = 1:numel(InterfacePoints)
    %        center = InterfacePoints(i);
    %        neighbors = getClosePoints(center, 0.004, gridHyPHM);
    %        normalSpeednew(neighbors,i) =normalSpeed(center);
    %     end
    %     normalSpeed = sum(normalSpeednew,2)./max(sum(normalSpeednew>0,2),1);
    imagesc(reshape(normalSpeed, 2 * 128 + 1, 128 + 1))
    %

    CFL = 1 / 4 * 1 / 10 * 1 / numPartitionsMicroscale / normalSpeedMax;
    %microscaleTimeStepSize = transportStepper.curtau / numMicroscaleTimeSteps;
    oldMicroscaleTime = currentTime;
    for j = 1:ceil(transportStepper.curtau/CFL)
        microscaleTimeStepSize = min(currentTime+transportStepper.curtau-oldMicroscaleTime, CFL);
        %             disp( ['    Level set substep ', num2str(j) ] );

        newMicroscaleTime = oldMicroscaleTime + microscaleTimeStepSize;
        % [] argument is unused in method (needed for implicit methods)
        currentLevelSetDataCells{1} = levelSetEquationTimeStep( ...
            newMicroscaleTime, ...
            oldMicroscaleTime, oldLevelSetDataCells{1}, microscaleGrid, normalSpeed, 1);
        oldMicroscaleTime = newMicroscaleTime;
        oldLevelSetDataCells{1} = currentLevelSetDataCells{1};

    end
    %   levelSet{1}(:,timeIterationStep+1) = reinitializeLevelSet( microscaleGrid, currentLevelSetDataCells{1});

    levelSet{1}(:, timeIterationStep + 1) = currentLevelSetDataCells{1};

    %     t_old = currentTime;
    %     currentTime = currentTime + dt_macro;
    %     oldLevelSetDataCells = currentLevelSetDataCells;
    currentTime = currentTime + macroscaleTimeStepSize;

    levelSetEvolutionTime(timeIterationStep) = toc;
    disp(['    ... done in ', ...
        num2str(levelSetEvolutionTime(timeIterationStep)), ' seconds.']);

    hydrogenTransport.L.setdata(timeIterationStep, levelSet{1}(:, timeIterationStep + 1));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Calculation of effective parameters

    disp('Calculation of effective parameters ...');
    tic;

    for i = 1:numberOfSlices

        [cellProblemSystemMatrix, rhs, isDoF, triangleVolumes, ...
            triangleSurfaces] = assembleCellProblem(microscaleGrid, ...
            levelSet{i}(:, timeIterationStep + 1));
        SOL = solveSystemFE(microscaleGrid, cellProblemSystemMatrix, rhs, isDoF);
        [diffusion, porosities] = computeDiffusionTensor(microscaleGrid, SOL, ...
            triangleVolumes);
        surfaceArea{i}(timeIterationStep + 1) = sum(triangleSurfaces);

    end
    Volume(timeIterationStep+1) = (0.05 * 0.1 - porosities);
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
    StokesL.computeLevel('s');

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


    %concentrations P_0(T)--> P_0(E)
    hydrogenDataE = zeros(gridHyPHM.numE, 1);


    for i = 1:gridHyPHM.numE
        if ((Levels(gridHyPHM.V0E(i, 1)) > -eps) & (Levels(gridHyPHM.V0E(i, 2)) > -eps))
            helper = gridHyPHM.T0E(i, :);
            hydrogenDataE(i, 1) = max(hydrogenData(helper));
        end
    end


    %%%%%%

    hydrogenTransportRhsData = zeros(gridHyPHM.numE, 1);
    hydrogenTransport.A.setdata(timeIterationStep, @(t, x) 1);
    hydrogenTransport.gF.setdata(timeIterationStep, hydrogenTransport.gF.getdata(0)+hydrogenTransportRhsData);
    hydrogenTransport.C.setdata(timeIterationStep, flow.getdata(1));


    SorceScale = zeros(gridHyPHM.numT, 1);
    for kT = 1:gridHyPHM.numT
        L = 0;
        for i = 1:3
            if Levels(gridHyPHM.V0E(gridHyPHM.E0T(kT, i), 1)) > -eps & Levels(gridHyPHM.V0E(gridHyPHM.E0T(kT, i), 2)) > -eps
                L = L + gridHyPHM.areaE(gridHyPHM.E0T(kT, i)) / gridHyPHM.areaT(kT);
            end
        end
        SorceScale(kT) = L;
    end

    speciesCells = {hydrogenTransport};
    nonlinearFunc = cell(1, 1);
    nonlinearFunc{1} = @(x, y) -SorceScale ...
        .* (x * rateCoefficientTST * 1000);
    % nonlinearFunc{2} = @(x,y) -nonlinearFunc{1}( x,y);

    nonlinearJacFunc = cell(1, 1);
    nonlinearJacFunc{1, 1} = @(x, y) -SorceScale ...
        .* (rateCoefficientTST * 1000);
    %         nonlinearJacFunc{1,2} = @(x,y) SorceScale * 0;
    %
    %         nonlinearJacFunc{2,1} = @(x,y) -nonlinearJacFunc{1,1}( x, y);
    %         nonlinearJacFunc{2,2} = @(x,y) nonlinearJacFunc{1,2}( x, y );
    newtonIteration(speciesCells, nonlinearFunc, nonlinearJacFunc, 2);
    %

    hydrogenTransport.U.setdata(timeIterationStep-1, hydrogenTransport.U.getdata(timeIterationStep - 1)-Corrector(:, 1));

    Rate(timeIterationStep+1) = -((hydrogenTransport.Q.getdata(timeIterationStep) - initialHydrogenConcentration * inletVelocity) .* (gridHyPHM.baryE(:, 1) > (0.1 - eps)))' * gridHyPHM.areaE / surfaceArea{1}(timeIterationStep);


end % while
% plot results related to figure 7
figure
plot(timeSteps(1:end), Rate);
xlabel('time [s]')
ylabel('average dissolution rate')

figure
plot(timeSteps(1:end), surfaceArea{1}(:));
xlabel('time [s]')
ylabel('surface Area')

figure
plot(timeSteps(1:end), Volume);
xlabel('time [s]')
ylabel('Grain Volume')
% hydrogenTransport.U.visualize;

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

function out = getInterfacePoints(Levels)

out = find((Levels < 0).*(Levels > -1 / (10 * 2 * 128)));

end

function out = getClosePoints(Center, dist, grid)
out = find(sqrt((grid.coordV(:, 1)-grid.coordV(Center, 1)).^2 + (grid.coordV(:, 2) - grid.coordV(Center, 2)).^2) < dist);
end
