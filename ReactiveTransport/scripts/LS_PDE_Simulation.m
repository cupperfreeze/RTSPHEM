%PDE level-set simulation
%cf. [0], Section 4.4 / 4.5

close all; %clear all;

%parameter
tic
global Solver
Solver = 'ilukIterative';
timer = zeros(3, 1); %Newton %Stokes %Level
%sparsityHit = [0.125*(1:12),1.6];
sparsityHit = [0.25 * (1:16)]; %time points to be hit (only ones saved)
sparsityHitCount = [0];
currentTime = 0;
reinitCount = 1;
finishTime = 4; %1.6;                                                        %final simulation time
numtimeSteps = 190; %90;                                                     %upper bound on number of time steps
saveConcentrationContent = zeros(numtimeSteps, 3); %Integrals of concentrations on fluid domain
FluxSave = zeros(numtimeSteps, 3); %Save flux values

%setup grid

%%
dimension = 2;
numPartitions = 199;
cellGrid = CartesianGrid(dimension, ...
    kron(ones(1, dimension), [-0.5, 0.5]), ...
    numPartitions*ones(1, dimension));
coord = cellGrid.coordinates;
temp = numPartitions + 1; %temporary variable for plotting
epsilon = 2 / numPartitions;
levelSave = zeros((numPartitions+1)^2, numtimeSteps, 'single'); %save level-set function at each step

timeSteps = linspace(currentTime, finishTime/(numtimeSteps - 1), 2);
timeStepsSave = linspace(currentTime, finishTime, numtimeSteps);
transportStepperSave = Stepper(timeStepsSave);
saveXi = zeros(temp*temp, numtimeSteps, 'uint8'); %save indicator function at each step
transportStepper = Stepper(timeSteps);
saveSolidArea = 0.2^2 * pi; %initial total solid surface area
saveDArea = saveSolidArea / 2; %initial surface area of D
gridHyPHMBasis = Grid(coord, cellGrid.triangles);

%Identify Edges
gridHyPHMBasis.idE(gridHyPHMBasis.baryE(:, 1) < -0.5+eps) = 4;
gridHyPHMBasis.idE(gridHyPHMBasis.baryE(:, 1) > 0.5-eps) = 2;
gridHyPHMBasis.idE(gridHyPHMBasis.baryE(:, 2) < -0.5+eps) = 1;
gridHyPHMBasis.idE(gridHyPHMBasis.baryE(:, 2) > 0.5-eps) = 3;

% Initial indicator function for different subdomains
Xi = ones(size(coord, 1), 1);
Xi((coord(:, 1) > 0) & (sqrt(coord(:, 1).^2 + coord(:, 2).^2) <= 0.2)) = 2;
Xi((coord(:, 1) <= 0) & (sqrt(coord(:, 1).^2 + coord(:, 2).^2) <= 0.2)) = 3;
saveXi(:, 1) = Xi;

%mineral densities
rhoP = 4;
rhoD = 20;

%Initial configuration Phi
initialPhiFunc = @(x) max(min(-(0.2 - epsilon - norm(x, 2)), 0.2 + epsilon - norm(x, 2)), (epsilon - norm(x(1), inf))*(abs(x(2)) <= 0.2)-100*(abs(x(2)) > 0.2));
coordCell = mat2cell(coord, ones(1, cellGrid.nodes), dimension);
Phi = cellfun(initialPhiFunc, coordCell);

%Calculate distance function d^i only up to this value
restrictDist = 3 * epsilon;

%Equilibrate Xi according to distFunctions
[distFunctions, Xi] = step1(Xi, Phi, cellGrid, restrictDist);
[dV, S] = step234(distFunctions, cellGrid, Xi, zeros(3), inf);

levelSet = dV;
levelSet(Xi == 1) = -levelSet(Xi == 1);
levelSet(Xi == 2) = max(levelSet(Xi == 2), 10^(-4));
levelSet(Xi == 3) = max(levelSet(Xi == 3), 10^(-4));


[gridHyPHM, levelSetNew] = localMeshRefinementTripel(gridHyPHMBasis, ... %obtain new grid aligned with current interface
    levelSet);

%Identify Edges
gridHyPHM.idE(gridHyPHM.baryE(:, 1) < -0.5+eps) = 4;
gridHyPHM.idE(gridHyPHM.baryE(:, 1) > 0.5-eps) = 2;
gridHyPHM.idE(gridHyPHM.baryE(:, 2) < -0.5+eps) = 1;
gridHyPHM.idE(gridHyPHM.baryE(:, 2) > 0.5-eps) = 3;


diffusionCoefficient = 0.2;
inletVelocity = 0; %inletVelocity=0 --> do not solve Stokes

saveConcentrationContent(1, 1) = 2 * (1 - 0.2^2 * pi);
saveConcentrationContent(1, 2) = 1 * (1 - 0.2^2 * pi);
saveConcentrationContent(1, 3) = 1 * (1 - 0.2^2 * pi);

% Define solute species and intial conditions
for i = 1:1
    Aconcentration = Variable(gridHyPHM, transportStepper, ...
        'A', 'P0');
    Aconcentration.setdata(2*ones(gridHyPHM.numT, 1));
    AconcentrationSave = Variable(gridHyPHM, transportStepperSave, ...
        'A', 'P0');
    ATransport.U = Aconcentration;


    BConcentration = Variable(gridHyPHM, transportStepper, ...
        'B', 'P0');
    BConcentration.setdata(1*ones(gridHyPHM.numT, 1));
    BconcentrationSave = Variable(gridHyPHM, transportStepperSave, ...
        'B', 'P0');

    BTransport.U = BConcentration;


    CConcentration = Variable(gridHyPHM, transportStepper, ...
        'C', 'P0');
    CConcentration.setdata(1*ones(gridHyPHM.numT, 1));
    CconcentrationSave = Variable(gridHyPHM, transportStepperSave, ...
        'C', 'P0');
    CTransport.U = CConcentration;


end


AconcentrationSave.setdata(0, ATransport.U.getdata(0));
BconcentrationSave.setdata(0, BTransport.U.getdata(0));
CconcentrationSave.setdata(0, CTransport.U.getdata(0));

%%

while transportStepperSave.next

    transportStepper.next;

    %Calculate distance function from each subdomain
    localTimer = tic;

    interfaceVelocities = [0, -1, -2; 0, 0, 0; 0, 0, 0];
    interfaceVelocities = interfaceVelocities - interfaceVelocities';

    BoundaryTriangles = abs(sum(levelSetNew(gridHyPHM.V0T(:, :)) > 0, 2)-2) < eps;
    indexBoundaryTriangles = find(BoundaryTriangles);
    baryBoundaryTriangles = gridHyPHM.baryT(BoundaryTriangles, :);
    localA = ATransport.U.getdata(1);
    localB = BTransport.U.getdata(1);
    localC = CTransport.U.getdata(1);

    %mark nodes at different fluid-solid interfaces
    [dVInitial, initialS] = step234Simple(distFunctions, cellGrid, Xi, interfaceVelocities);
    IdInterface1 = abs(abs(initialS)-1) < eps;
    IdInterface2 = abs(abs(initialS)-2) < eps;

    interfaceVelocitiesNodes = zeros(gridHyPHMBasis.numV, 1);

    for i = find(~isinf(initialS))'
        [~, idx] = min(2*abs(baryBoundaryTriangles(:, 1) - gridHyPHMBasis.coordV(i, 1))+ ...
            abs (baryBoundaryTriangles(:, 2) - gridHyPHMBasis.coordV(i, 2)));
        %    [~,idx] =min( sqrt((baryBoundaryTriangles(:,1) - gridHyPHMBasis.coordV(i,1)).^2 + ...
        %                        + (baryBoundaryTriangles(:,2) - gridHyPHMBasis.coordV(i,2)).^2));
        interfaceVelocitiesNodes(i) = -(localB(indexBoundaryTriangles(idx)) ./ localA(indexBoundaryTriangles(idx)) - 1) .* IdInterface1(i) ...
            -(localB(indexBoundaryTriangles(idx)) .* localC(indexBoundaryTriangles(idx)) - 1) .* IdInterface2(i);
    end


    saveXi(:, transportStepperSave.curstep+1) = Xi;
    initialS(~isinf(dVInitial)) = -interfaceVelocitiesNodes(~isinf(dVInitial)) .* sign(initialS(~isinf(dVInitial)));


    %Reconstruct Vonoroi Interface and velocity extension
    [dV, S] = reinitializeLevelSetWithS(cellGrid, dVInitial, 0.2, initialS');
    CFL = min(0.4/max(S(~isinf(S)))/numPartitions, 0.04);


    %update Time stepping accoding to CFL
    transportStepperSave.setTimeStepSize(transportStepperSave.curstep, min(CFL, sparsityHit(numel(sparsityHitCount)) - currentTime));
    transportStepper.setTimeStepSize(transportStepper.curstep, min(CFL, sparsityHit(numel(sparsityHitCount)) - currentTime));


    %evolve geometry
    for j = 1:ceil(transportStepperSave.curtau/CFL)
        StepSizeTime = min(CFL, transportStepperSave.curtime-currentTime);
        currentTime = currentTime + StepSizeTime;
        if abs(sparsityHit(numel(sparsityHitCount))-currentTime) < 10 * eps
            sparsityHitCount = [sparsityHitCount, transportStepperSave.curstep];
        end

        %Evolve Phi
        Phi = levelSetEquationTimeStep(StepSizeTime, 0, Phi, ...
            cellGrid, S', 1);

        reinitCount = reinitCount + 1;
    end

    %Perform reinitialization step
    if mod(reinitCount, 30) == 0
        Phi = epsilon - dV;
    end

    %Updated Xi according to Phi
    [distFunctions, Xi] = step1(Xi, Phi, cellGrid, restrictDist);
    [dV, ~] = step234(distFunctions, cellGrid, Xi, zeros(3), 0.2);

    % prepare levelSet such that {levelSet=0} approximated total
    % fluid-solid interface
    levelSet = dV;
    levelSet(Xi == 1) = -levelSet(Xi == 1);
    levelSet(Xi == 2) = max(levelSet(Xi == 2), 10^(-4));
    levelSet(Xi == 3) = max(levelSet(Xi == 3), 10^(-4));

    % Evaluate surface areas
    [~, ~, ~, triangleVolumes, ~] ...
        = assembleCellProblem(cellGrid, levelSet);
    saveSolidArea = [saveSolidArea; 1 - sum(triangleVolumes)];

    levelTemp = dV;
    levelTemp(Xi ~= 2) = -levelTemp(Xi ~= 2);
    [~, ~, ~, triangleVolumes, ~] ...
        = assembleCellProblem(cellGrid, levelTemp);
    saveDArea = [saveDArea; 1 - sum(triangleVolumes)];
    [gridHyPHM, levelSetNew] = localMeshRefinementTripel(gridHyPHMBasis, ...
        levelSet);

    levelSave(:, reinitCount) = levelSet;

    %evaluate length of interfaces and coordinates of triple points
    [interfaceLength{reinitCount}, coordTriple{reinitCount}] = evaluateInterface(cellGrid, Xi, distFunctions, true);
    timer(3) = timer(3) + toc(localTimer);

    %Identify Edges
    gridHyPHM.idE(gridHyPHM.baryE(:, 1) < -0.5+eps) = 4;
    gridHyPHM.idE(gridHyPHM.baryE(:, 1) > 0.5-eps) = 2;
    gridHyPHM.idE(gridHyPHM.baryE(:, 2) < -0.5+eps) = 1;
    gridHyPHM.idE(gridHyPHM.baryE(:, 2) > 0.5-eps) = 3;

    %% Don't allow for inclusions

    aux = gridHyPHM.V0E(:, :);
    numNeighbors = zeros(gridHyPHM.numV, 1);
    numSolidNeighbors = zeros(gridHyPHM.numV, 1);
    for i = 1:size(aux, 1)
        numNeighbors(aux(i, 1)) = numNeighbors(aux(i, 1)) + 1;
        numNeighbors(aux(i, 2)) = numNeighbors(aux(i, 2)) + 1;
        numSolidNeighbors(aux(i, 1)) = numSolidNeighbors(aux(i, 1)) + (levelSetNew(aux(i, 2)) > 0);
        numSolidNeighbors(aux(i, 2)) = numSolidNeighbors(aux(i, 2)) + (levelSetNew(aux(i, 1)) > 0);
    end


    levelSetNew(find(levelSetNew < 0 & (numSolidNeighbors == numNeighbors) & (numNeighbors) > 0)) = eps;

    %% Setup fluid simulation
    localTimer = tic;
    phead = Variable(gridHyPHM, transportStepper, 'pressure head', 'P1');
    StokesVelocity = Variable(gridHyPHM, transportStepper, 'Stokes Velocity', 'P2P2');
    phead.setdata(0, zeros(gridHyPHM.numV, 1));
    StokesVelocity.setdata(zeros(gridHyPHM.numV + gridHyPHM.numE, 2));

    StokesL = StokesLEVEL(gridHyPHM, transportStepper, 'Stokes problem');
    StokesL.id2D = {4, 3, 1};
    StokesL.id2N = {2};
    % StokesL.uD.setdata(@(t, x) inletVelocity*(x(1)<-0.5+eps)*[1;0] );
    % StokesL.uD.setdata([inletVelocity, 0] .* ( [gridHyPHM.coordV(:,1),gridHyPHM.coordV(:,1);gridHyPHM.baryE(:,1),gridHyPHM.baryE(:,1)] < -0.5+eps ));
    StokesL.uD.setdata([inletVelocity, 0].*([gridHyPHM.coordV(:, 1), gridHyPHM.coordV(:, 1); gridHyPHM.baryE(:, 1), gridHyPHM.baryE(:, 1)] < -0.5 + eps).* ...
        [1 - (2 * gridHyPHM.coordV(:, 2)).^2; 1 - (2 * gridHyPHM.baryE(:, 2)).^2]);
    StokesL.F.setdata(zeros(gridHyPHM.numV + gridHyPHM.numE, 2));
    StokesL.U = StokesVelocity;
    StokesL.P = phead;

    StokesL.L.setdata(levelSetNew);
    if abs(inletVelocity) > 0 %only solve when necessary
        StokesL.computeLevel('s');
    end
    timer(2) = timer(2) + toc(localTimer);

    % Setup transport problem on adapted grid, specify exterior boundary conditions
    for i = 1:1
        temp = StokesL.U.getdata(1);
        flow = temp(gridHyPHM.numV+1:end, 1) .* gridHyPHM.nuE(:, 1) + temp(gridHyPHM.numV+1:end, 2) .* gridHyPHM.nuE(:, 2);

        Aconcentration = Variable(gridHyPHM, transportStepper, ...
            'A', 'P0');

        ATransport = TransportLEVEL(gridHyPHM, transportStepper, 'A Transport');
        ATransport.D.setdata(diffusionCoefficient*eye(2));
        ATransport.id2N = {1, 2, 3};
        ATransport.id2F = {4};
        ATransport.U = Aconcentration;
        ATransport.gF = Variable(gridHyPHM, transportStepper, 'gF', 'P0E');
        %ATransport.gF.setdata(  - 2*inletVelocity * ( gridHyPHM.baryE(:,1) < -0.5+eps ) );
        %ATransport.gF.setdata(  - 2*inletVelocity * ( gridHyPHM.baryE(:,1) < -0.5+eps ).* (1-(2*gridHyPHM.baryE(:,2)).^2) );
        ATransport.gF.setdata(-1.4056*inletVelocity*(gridHyPHM.baryE(:, 1) < -0.5 + eps).*(1 - (2 * gridHyPHM.baryE(:, 2)).^2));
        ATransport.A = Variable(gridHyPHM, transportStepper, 'Porosity', 'P0');
        ATransport.A.setdata(ones(gridHyPHM.numT, 1));
        ATransport.C.setdata(flow);
        ATransport.isUpwind = 'exp';


        BConcentration = Variable(gridHyPHM, transportStepper, ...
            'B', 'P0');


        BTransport = TransportLEVEL(gridHyPHM, transportStepper, 'B Transport');
        BTransport.D.setdata(diffusionCoefficient*eye(2));
        BTransport.id2N = {1, 2, 3};
        BTransport.id2F = {4};
        BTransport.U = BConcentration;
        BTransport.gF = Variable(gridHyPHM, transportStepper, 'gF', 'P0E');
        %BTransport.gF.setdata( - 1*inletVelocity * ( gridHyPHM.baryE(:,1) < -0.5+eps ) );
        %BTransport.gF.setdata( - 1*inletVelocity *  ( gridHyPHM.baryE(:,1) < -0.5+eps ).* (1-(2*gridHyPHM.baryE(:,2)).^2));
        BTransport.gF.setdata(-1.4056*inletVelocity*(gridHyPHM.baryE(:, 1) < -0.5 + eps).*(1 - (2 * gridHyPHM.baryE(:, 2)).^2));
        BTransport.A = Variable(gridHyPHM, transportStepper, 'Porosity', 'P0');
        BTransport.A.setdata(ones(gridHyPHM.numT, 1));
        BTransport.C.setdata(flow);
        BTransport.isUpwind = 'exp';

        CConcentration = Variable(gridHyPHM, transportStepper, ...
            'C', 'P0');

        CTransport = TransportLEVEL(gridHyPHM, transportStepper, 'C Transport');
        CTransport.D.setdata(diffusionCoefficient*eye(2));
        CTransport.id2N = {1, 2, 3};
        CTransport.id2F = {4};
        CTransport.U = CConcentration;
        CTransport.gF = Variable(gridHyPHM, transportStepper, 'gF', 'P0E');
        %CTransport.gF.setdata( - 1*inletVelocity * ( gridHyPHM.baryE(:,1) < -0.5+eps ) );
        %CTransport.gF.setdata( - 1*inletVelocity * ( gridHyPHM.baryE(:,1) < -0.5+eps ).* (1-(2*gridHyPHM.baryE(:,2)).^2) );
        CTransport.gF.setdata(-0.7114*inletVelocity*(gridHyPHM.baryE(:, 1) < -0.5 + eps).*(1 - (2 * gridHyPHM.baryE(:, 2)).^2));
        CTransport.A = Variable(gridHyPHM, transportStepper, 'Porosity', 'P0');
        CTransport.A.setdata(ones(gridHyPHM.numT, 1));
        CTransport.C.setdata(flow);
        CTransport.isUpwind = 'exp';
    end


    %interpolate data from old to new grid
    now = transportStepperSave.curstep - 1;
    disp(['Step: ', num2str(now)]);
    temp = AconcentrationSave.getTSI(now);
    ATransport.U.setdata(0, temp(gridHyPHM.baryT(:, 1), gridHyPHM.baryT(:, 2)));
    temp = BconcentrationSave.getTSI(now);
    BTransport.U.setdata(0, temp(gridHyPHM.baryT(:, 1), gridHyPHM.baryT(:, 2)));
    temp = CconcentrationSave.getTSI(now);
    CTransport.U.setdata(0, temp(gridHyPHM.baryT(:, 1), gridHyPHM.baryT(:, 2)));

    ATransport.L.setdata(levelSetNew);
    BTransport.L.setdata(levelSetNew);
    CTransport.L.setdata(levelSetNew);

    %Determine triangles at interior boundary and strength of source term
    SorceScale = zeros(gridHyPHM.numT, 1);
    IdInterface1 = zeros(gridHyPHMBasis.numT, 1);
    IdInterface2 = zeros(gridHyPHMBasis.numT, 1);

    IdInterface1(any(Xi(gridHyPHMBasis.V0T(:, :)) == 1, 2) & any(Xi(gridHyPHMBasis.V0T(:, :)) == 2, 2)) = 1;
    IdInterface2(any(Xi(gridHyPHMBasis.V0T(:, :)) == 1, 2) & any(Xi(gridHyPHMBasis.V0T(:, :)) == 3, 2)) = 1;

    temp = scatteredInterpolant(gridHyPHMBasis.baryT(:, 1), gridHyPHMBasis.baryT(:, 2), IdInterface1, 'linear');
    IdInterface1 = ceil(temp(gridHyPHM.baryT(:, 1), gridHyPHM.baryT(:, 2)));
    temp = scatteredInterpolant(gridHyPHMBasis.baryT(:, 1), gridHyPHMBasis.baryT(:, 2), IdInterface2, 'linear');
    IdInterface2 = ceil(temp(gridHyPHM.baryT(:, 1), gridHyPHM.baryT(:, 2)));


    for kT = find(IdInterface1 | IdInterface2)'
        L = 0;
        for i = 1:3
            if levelSetNew(gridHyPHM.V0E(gridHyPHM.E0T(kT, i), 1)) > -eps & levelSetNew(gridHyPHM.V0E(gridHyPHM.E0T(kT, i), 2)) > -eps
                L = L + gridHyPHM.areaE(gridHyPHM.E0T(kT, i)) / gridHyPHM.areaT(kT);
            end
        end
        SorceScale(kT) = L;
    end

    %Define sorce term (reaction rates)
    nonlinearFunc = cell(3, 1);
    nonlinearFunc{1} = @(a, b, c) (IdInterface1 .* (b ./ max(a, sqrt(eps)) - 1) .* (rhoD + a) + IdInterface2 .* (b .* c - 1) .* a) .* SorceScale;
    nonlinearFunc{2} = @(a, b, c) (-IdInterface1 .* (b ./ max(a, sqrt(eps)) - 1) .* (rhoD - b) - IdInterface2 .* (b .* c - 1) .* (rhoP - b)) .* SorceScale;
    nonlinearFunc{3} = @(a, b, c) (IdInterface1 .* (b ./ max(a, sqrt(eps)) - 1) .* c - IdInterface2 .* (b .* c - 1) .* (rhoP - c)) .* SorceScale;

    %Analytical derivatives for Jacobian
    nonlinearJacFunc = cell(3, 3);
    nonlinearJacFunc{1, 1} = @(a, b, c) (-IdInterface1 .* (b ./ max(a.^2, sqrt(eps))) .* (rhoD + a) + IdInterface1 .* (b ./ max(a, sqrt(eps)) - 1) ...
        +IdInterface2 .* (b .* c - 1)) .* SorceScale;
    nonlinearJacFunc{1, 2} = @(a, b, c) (IdInterface1 .* (1 ./ max(a, sqrt(eps))) .* (rhoD + a) + IdInterface2 .* c .* a) .* SorceScale;
    nonlinearJacFunc{1, 3} = @(a, b, c) IdInterface2 .* b .* a .* SorceScale;

    nonlinearJacFunc{2, 1} = @(a, b, c) (IdInterface1 .* (b ./ max(a.^2, sqrt(eps))) .* (rhoD - b)) .* SorceScale;
    nonlinearJacFunc{2, 2} = @(a, b, c) (-IdInterface2 .* c .* (rhoP - b) + IdInterface2 .* (b .* c - 1) ...
        -IdInterface1 .* (1 ./ max(a, sqrt(eps))) .* (rhoD - b) + IdInterface1 .* (b ./ max(a, sqrt(eps)))) .* SorceScale;
    nonlinearJacFunc{2, 3} = @(a, b, c) -IdInterface2 .* b .* (rhoP - b) .* SorceScale;

    nonlinearJacFunc{3, 1} = @(a, b, c) -IdInterface1 .* (b ./ max(a.^2, sqrt(eps))) .* c .* SorceScale;
    nonlinearJacFunc{3, 2} = @(a, b, c) (IdInterface1 .* 1 ./ max(a, sqrt(eps)) .* c - IdInterface2 .* c .* (rhoP - c)) .* SorceScale;
    nonlinearJacFunc{3, 3} = @(a, b, c) (IdInterface1 .* (b ./ max(a, sqrt(eps)) - 1) - IdInterface2 .* b .* (rhoP - c) + IdInterface2 .* (b .* c - 1)) .* SorceScale;

    %Call non-linear solver
    speciesCells = {ATransport; BTransport; CTransport};
    localTimer = tic;
    newtonIteration(speciesCells, nonlinearFunc, nonlinearJacFunc, 20);
    timer(1) = timer(1) + toc(localTimer);

    %save concentrations and grid of current step
    AconcentrationSave.setdataGrid(transportStepperSave.curstep, ATransport.U.getdata(1), gridHyPHM);
    BconcentrationSave.setdataGrid(transportStepperSave.curstep, BTransport.U.getdata(1), gridHyPHM);
    CconcentrationSave.setdataGrid(transportStepperSave.curstep, CTransport.U.getdata(1), gridHyPHM);

    % to save memory, delete unnecessary time steps
    if sparsityHitCount(end) == transportStepperSave.curstep
        for i = (sparsityHitCount(end -1) + 1):(sparsityHitCount(end) - 1)
            AconcentrationSave.grids{i+1} = 0;
            k = AconcentrationSave.grid.numT;
            AconcentrationSave.setdata(i, sparse(k, 1));
            BconcentrationSave.grids{i+1} = 0;
            k = BconcentrationSave.grid.numT;
            BconcentrationSave.setdata(i, sparse(k, 1));
            CconcentrationSave.grids{i+1} = 0;
            k = CconcentrationSave.grid.numT;
            CconcentrationSave.setdata(i, sparse(k, 1));
        end
    end

    % save flux variable at inflow and outflow boundary
    if transportStepperSave.curstep == 1
        FluxSave(1, 1) = (sum(ATransport.Q.getdata(1) .* gridHyPHM.areaE .* (abs(gridHyPHM.baryE(:, 1)) > 0.5 - eps)));
        FluxSave(1, 2) = (sum(BTransport.Q.getdata(1) .* gridHyPHM.areaE .* (abs(gridHyPHM.baryE(:, 1)) > 0.5 - eps)));
        FluxSave(1, 3) = (sum(CTransport.Q.getdata(1) .* gridHyPHM.areaE .* (abs(gridHyPHM.baryE(:, 1)) > 0.5 - eps)));
    end

    FluxSave(transportStepperSave.curstep+1, 1) = single(sum(ATransport.Q.getdata(1) .* gridHyPHM.areaE .* (abs(gridHyPHM.baryE(:, 1)) > 0.5 - eps)));
    FluxSave(transportStepperSave.curstep+1, 2) = single(sum(BTransport.Q.getdata(1) .* gridHyPHM.areaE .* (abs(gridHyPHM.baryE(:, 1)) > 0.5 - eps)));
    FluxSave(transportStepperSave.curstep+1, 3) = single(sum(CTransport.Q.getdata(1) .* gridHyPHM.areaE .* (abs(gridHyPHM.baryE(:, 1)) > 0.5 - eps)));

    saveConcentrationContent(transportStepperSave.curstep+1, 1) = sum(ATransport.U.getdata(1).*gridHyPHM.areaT);
    saveConcentrationContent(transportStepperSave.curstep+1, 2) = sum(BTransport.U.getdata(1).*gridHyPHM.areaT);
    saveConcentrationContent(transportStepperSave.curstep+1, 3) = sum(CTransport.U.getdata(1).*gridHyPHM.areaT);


    transportStepper.prev;
    if (currentTime >= finishTime - eps) | transportStepperSave.curstep == numtimeSteps
        break
    end

end
toc
%plot evolved functions
temp = numPartitions + 1;
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

%[interfaceLength, coordTriple] = evaluateInterface(cellGrid, Xi, distFunctions, false)