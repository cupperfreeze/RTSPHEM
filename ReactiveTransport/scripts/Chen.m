% MicroScale simulation using Voronoi Implicit Interface Method (VIIM).
% Setup inspired by paper [Chen, 2014].
% cf. [2] Section 4.2
close all
%parameter

global Solver
Solver = 'ilukIterative';
currentTime = 0;
reinitCount = 0;
finishTime = 15;
numtimeSteps = 200; %upper bound on timeSteps

%setup grid

%%
dimension = 2;
lengthX = 0.14;
lengthY = 0.09;
numPartitions = 210;
cellGrid = CartesianGrid(dimension, ...
    [0, lengthX, 0, lengthY], ...
    numPartitions*[1, 90 / 140]);
coord = cellGrid.coordinates;
temp = numPartitions + 1; %temporary variable for plotting
epsilon = sqrt(2) * lengthX / numPartitions; %epsilon in VIIM
restrictDist = 3 * epsilon; %cut-off distance in VIIM

% parameter
k1 = 10^(-7);
k2 = 10^(-9);
K1 = 10^(-2);
K2 = 1.25 * 10^(13);
molarVolume = 3 * 10^4;
diffusionCoefficient = 10^(-3);
inletVelocity = 0;

timeSteps = linspace(currentTime, finishTime/(numtimeSteps - 1), 2);
timeStepsSave = linspace(currentTime, finishTime, numtimeSteps);
transportStepperSave = Stepper(timeStepsSave);
saveXi = zeros((numPartitions+1)*(round(numPartitions * 9 / 14 + 1)), numtimeSteps);
transportStepper = Stepper(timeSteps);
gridHyPHMBasis = Grid(coord, cellGrid.triangles);
gridHyPHMBasis = localMeshRefinementTripel(gridHyPHMBasis, ones(gridHyPHMBasis.numT, 1));

%Identify Edges
gridHyPHMBasis.idE(gridHyPHMBasis.baryE(:, 1) < eps) = 4;
gridHyPHMBasis.idE(gridHyPHMBasis.baryE(:, 1) > lengthX-eps) = 2;
gridHyPHMBasis.idE(gridHyPHMBasis.baryE(:, 2) < eps) = 1;
gridHyPHMBasis.idE(gridHyPHMBasis.baryE(:, 2) > lengthY-eps) = 3;

% Initial indicator function for different subdomains
Xi = ones(size(coord, 1), 1);
Xi(((coord(:, 2) > 0.060) | (coord(:, 2) < 0.030)) & 0.01 < coord(:, 1) & coord(:, 1) < 0.13) = 2;
X = [repmat(linspace(0.02, 0.12, 10), 1, 6)];
Y = [repmat(0.00, 1, 10), repmat(0.01, 1, 10), repmat(0.02, 1, 10), repmat(0.07, 1, 10), repmat(0.08, 1, 10), repmat(0.09, 1, 10)];
radius = 0.003;
% for i=1:numel(X)
% Xi(sqrt((coord(:,1)-X(i)).^2 +(coord(:,2)-Y(i)).^2)<radius) = 3;
% end
%saveXi(:,1) = Xi;


%Initial configuration Phi
initialPhiFunc = @(x) min(max(.060 + epsilon - norm([x(1) - 0.070, x(2) - .120], inf), .060 + epsilon - norm([x(1) - .070, x(2) + .030], inf)), ...
    min(-.060 + epsilon + norm([x(1) - .070, x(2) - .120], inf), -.060 + epsilon + norm([x(1) - .070, x(2) + .030], inf)));
% initialPhiFunc = drawcircle(X,Y,radius, epsilon,initialPhiFunc);
[initialPhiFunc, Xi] = drawcircle(X, Y, repmat(radius, 1, 60), epsilon, initialPhiFunc, Xi, coord);
saveXi(:, 1) = Xi;
coordCell = mat2cell(coord, ones(1, cellGrid.nodes), dimension);
Phi = cellfun(initialPhiFunc, coordCell);
surf(reshape(Phi, numPartitions + 1, round(numPartitions * 9 / 14 + 1)));
%Calculate distance function d^i only up to this value


%Equilibrate Xi according to distFunctions
[distFunctions, Xi] = step1(Xi, Phi, cellGrid, restrictDist);
[dV, S] = step234(distFunctions, cellGrid, Xi, zeros(3), restrictDist);

levelSet = dV;
levelSet(Xi == 1) = -levelSet(Xi == 1);
%     levelSet(Xi==2) = max(levelSet(Xi==2),10^(-4));
%     levelSet(Xi==3) = max(levelSet(Xi==3),10^(-4));


[gridHyPHM, levelSetNew] = localMeshRefinementTripel(gridHyPHMBasis, ...
    levelSet);

%Identify Edges
gridHyPHM.idE(gridHyPHM.baryE(:, 1) < eps) = 4;
gridHyPHM.idE(gridHyPHM.baryE(:, 1) > lengthX-eps) = 2;
gridHyPHM.idE(gridHyPHM.baryE(:, 2) < eps) = 1;
gridHyPHM.idE(gridHyPHM.baryE(:, 2) > lengthY-eps) = 3;

% Initial conditions
for i = 1:1
    Aconcentration = Variable(gridHyPHM, transportStepper, ...
        'A', 'P0');
    Aconcentration.setdata(2*10^(-6)*ones(gridHyPHM.numT, 1));
    AconcentrationSave = Variable(gridHyPHM, transportStepperSave, ...
        'A', 'P0');
    ATransport.U = Aconcentration;


    BConcentration = Variable(gridHyPHM, transportStepper, ...
        'B', 'P0');
    BConcentration.setdata(2*10^(-4)*ones(gridHyPHM.numT, 1));
    BconcentrationSave = Variable(gridHyPHM, transportStepperSave, ...
        'B', 'P0');

    BTransport.U = BConcentration;


    CConcentration = Variable(gridHyPHM, transportStepper, ...
        'C', 'P0');
    CConcentration.setdata(zeros(gridHyPHM.numT, 1));
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

    %Identify different fluid-solid interface traingles

    interfaceVelocities = [0, -1, -2; 0, 0, 0; 0, 0, 0];
    interfaceVelocities = interfaceVelocities - interfaceVelocities';

    BoundaryTriangles = abs(sum(levelSetNew(gridHyPHM.V0T(:, :)) > -eps, 2)-2) < eps;
    indexBoundaryTriangles = find(BoundaryTriangles);
    baryBoundaryTriangles = gridHyPHM.baryT(BoundaryTriangles, :);
    localA = ATransport.U.getdata(1);
    localB = BTransport.U.getdata(1);
    localC = CTransport.U.getdata(1);


    [dVInitial, initialS] = step234Simple(distFunctions, cellGrid, Xi, interfaceVelocities);
    IdInterface1 = abs(abs(initialS)-1) < eps;
    IdInterface2 = abs(abs(initialS)-2) < eps;

    interfaceVelocitiesNodes = zeros(gridHyPHMBasis.numV, 1);

    %calculate normal interface velocity
    for i = find(~isinf(initialS))'
        [~, idx] = min(2*abs(baryBoundaryTriangles(:, 1) - gridHyPHMBasis.coordV(i, 1))+ ...
            abs (baryBoundaryTriangles(:, 2) - gridHyPHMBasis.coordV(i, 2)));
        interfaceVelocitiesNodes(i) = -molarVolume * k1 * (K1 * localB(indexBoundaryTriangles(idx)) ./ localA(indexBoundaryTriangles(idx)) - 1) .* IdInterface1(i) ...
            -0.25 * molarVolume * k2 * (K2 * localB(indexBoundaryTriangles(idx)) .* localC(indexBoundaryTriangles(idx)) - 1) .* IdInterface2(i);
    end


    saveXi(:, transportStepperSave.curstep+1) = Xi;
    initialS(~isinf(dVInitial)) = -interfaceVelocitiesNodes(~isinf(dVInitial)) .* sign(initialS(~isinf(dVInitial)));


    %Reconstruct Vonoroi Interface and velocity extension
    [dV, S] = reinitializeLevelSetWithS(cellGrid, dVInitial, inf, initialS');
    maxS = max(S(~isinf(S)));
    CFL = abs(lengthX*0.5/maxS/numPartitions);
    if transportStepperSave.curstep < 5
        CFL = min(CFL, 0.1);
    end
    %update Time stepping accoding to CFL
    transportStepperSave.setTimeStepSize(transportStepperSave.curstep, min(CFL, finishTime - currentTime));
    transportStepper.setTimeStepSize(transportStepper.curstep, min(CFL, finishTime - currentTime));

    for j = 1:ceil(transportStepperSave.curtau/CFL)
        StepSizeTime = min(CFL, transportStepperSave.curtime-currentTime);
        currentTime = currentTime + StepSizeTime;

        %Evolve Phi
        Phi = levelSetEquationTimeStep(StepSizeTime, 0, Phi, ...
            cellGrid, S');

        reinitCount = reinitCount + 1;
    end

    if mod(reinitCount, 30) == 0
        Phi = epsilon - dV;
    end

    %Updated Xi according to Phi
    [distFunctions, Xi] = step1(Xi, Phi, cellGrid, restrictDist);
    [dV, ~] = step234(distFunctions, cellGrid, Xi, zeros(3), restrictDist);

    levelSet = dV;
    levelSet(Xi == 1) = -levelSet(Xi == 1);
    levelSet(Xi == 2) = max(levelSet(Xi == 2), 10^(-4));
    levelSet(Xi == 3) = max(levelSet(Xi == 3), 10^(-4));


    [gridHyPHM, levelSetNew] = localMeshRefinementTripel(gridHyPHMBasis, ...
        levelSet);

    %Identify Edges
    gridHyPHM.idE(gridHyPHM.baryE(:, 1) < eps) = 4;
    gridHyPHM.idE(gridHyPHM.baryE(:, 1) > lengthX-eps) = 2;
    gridHyPHM.idE(gridHyPHM.baryE(:, 2) < eps) = 1;
    gridHyPHM.idE(gridHyPHM.baryE(:, 2) > lengthY-eps) = 3;

    %% Don't allow for singular inclusions

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


    % calculate Stokes velocity and pressure field if necessary
    phead = Variable(gridHyPHM, transportStepper, 'pressure head', 'P1');
    StokesVelocity = Variable(gridHyPHM, transportStepper, 'Stokes Velocity', 'P2P2');
    phead.setdata(0, zeros(gridHyPHM.numV, 1));
    StokesVelocity.setdata(zeros(gridHyPHM.numV + gridHyPHM.numE, 2));

    StokesL = StokesLEVEL(gridHyPHM, transportStepper, 'Stokes problem');
    StokesL.id2D = {4, 3, 1};
    StokesL.id2N = {2};
    StokesL.uD.setdata(@(t, x) inletVelocity*(x(1) < eps)*[1; 0]);
    StokesL.F.setdata(zeros(gridHyPHM.numV + gridHyPHM.numE, 2));
    StokesL.U = StokesVelocity;
    StokesL.P = phead;

    StokesL.L.setdata(levelSetNew);
    if abs(inletVelocity) > 0
        StokesL.computeLevel('s');
    end

    % setup new variables on adapted mesh, set boundary conditions
    for i = 1:1
        temp = StokesL.U.getdata(1);
        flow = temp(gridHyPHM.numV+1:end, 1) .* gridHyPHM.nuE(:, 1) + temp(gridHyPHM.numV+1:end, 2) .* gridHyPHM.nuE(:, 2);

        Aconcentration = Variable(gridHyPHM, transportStepper, ...
            'A', 'P0');

        ATransport = TransportLEVEL(gridHyPHM, transportStepper, 'A Transport');
        ATransport.id2D = {4};
        ATransport.uD.setdata(2*10^(-6)*ones(gridHyPHM.numE, 1));
        ATransport.D.setdata(diffusionCoefficient*eye(2));
        ATransport.id2N = {1, 3};
        ATransport.id2F = {2};
        ATransport.U = Aconcentration;
        % calciumTransport.uD = Variable( g_HyPHM, st, 'uD', 'P0E' );
        % calciumTransport.uD.setdata( zeros( g_HyPHM.numE, 1 ) );
        ATransport.gF = Variable(gridHyPHM, transportStepper, 'gF', 'P0E');
        ATransport.gF.setdata(-2*inletVelocity*(gridHyPHM.baryE(:, 1) < eps));
        ATransport.A = Variable(gridHyPHM, transportStepper, 'Porosity', 'P0');
        ATransport.A.setdata(ones(gridHyPHM.numT, 1));
        ATransport.C.setdata(flow);
        ATransport.isUpwind = 'exp';


        BConcentration = Variable(gridHyPHM, transportStepper, ...
            'B', 'P0');


        % carbonateConcentration.setdata( 0, @(t,x) sqrt( equilibriumRateConstant ) );
        BTransport = TransportLEVEL(gridHyPHM, transportStepper, 'B Transport');
        BTransport.id2D = {4};
        BTransport.uD.setdata(zeros(gridHyPHM.numE, 1));
        BTransport.D.setdata(diffusionCoefficient*eye(2));
        BTransport.id2N = {1, 3};
        BTransport.id2F = {2};
        BTransport.U = BConcentration;
        % carbonateTransport.uD = Variable( g_HyPHM, st, 'uD', 'P0E' );
        % carbonateTransport.uD.setdata( zeros( g_HyPHM.numE, 1 ) );
        BTransport.gF = Variable(gridHyPHM, transportStepper, 'gF', 'P0E');
        BTransport.gF.setdata(-1*inletVelocity*(gridHyPHM.baryE(:, 1) < eps));
        BTransport.A = Variable(gridHyPHM, transportStepper, 'Porosity', 'P0');
        BTransport.A.setdata(ones(gridHyPHM.numT, 1));
        BTransport.C.setdata(flow);
        BTransport.isUpwind = 'exp';

        CConcentration = Variable(gridHyPHM, transportStepper, ...
            'C', 'P0');

        % magnesiumConcentration.setdata( 0, @(t,x) sqrt( equilibriumRateConstant ) );
        CTransport = TransportLEVEL(gridHyPHM, transportStepper, 'C Transport');
        CTransport.id2D = {4};
        CTransport.uD.setdata(2*10^(-6)*ones(gridHyPHM.numE, 1));
        CTransport.D.setdata(diffusionCoefficient*eye(2));
        CTransport.id2N = {1, 3};
        CTransport.id2F = {2};
        CTransport.U = CConcentration;
        % magnesiumTransport.uD = Variable( g_HyPHM, st, 'uD', 'P0E' );
        % magnesiumTransport.uD.setdata( zeros( g_HyPHM.numE, 1 ) );
        CTransport.gF = Variable(gridHyPHM, transportStepper, 'gF', 'P0E');
        CTransport.gF.setdata(-1*inletVelocity*(gridHyPHM.baryE(:, 1) < eps));
        CTransport.A = Variable(gridHyPHM, transportStepper, 'Porosity', 'P0');
        CTransport.A.setdata(ones(gridHyPHM.numT, 1));
        CTransport.C.setdata(flow);
        CTransport.isUpwind = 'exp';
    end %Initialize new Transport instances


    %set data on new mesh
    now = transportStepperSave.curstep - 1;
    temp = AconcentrationSave.getTSI(now);
    ATransport.U.setdata(0, temp(gridHyPHM.baryT(:, 1), gridHyPHM.baryT(:, 2)));
    temp = BconcentrationSave.getTSI(now);
    BTransport.U.setdata(0, temp(gridHyPHM.baryT(:, 1), gridHyPHM.baryT(:, 2)));
    temp = CconcentrationSave.getTSI(now);
    CTransport.U.setdata(0, temp(gridHyPHM.baryT(:, 1), gridHyPHM.baryT(:, 2)));

    %restrict to fluid domain
    ATransport.L.setdata(levelSetNew);
    BTransport.L.setdata(levelSetNew);
    CTransport.L.setdata(levelSetNew);

    %identify interfaces and local source strength
    SorceScale = zeros(gridHyPHM.numT, 1);
    IdInterface1 = zeros(gridHyPHMBasis.numT, 1);
    IdInterface2 = zeros(gridHyPHMBasis.numT, 1);

    IdInterface1(any(Xi(gridHyPHMBasis.V0T(:, :)) == 1, 2) & any(Xi(gridHyPHMBasis.V0T(:, :)) == 2, 2)) = 1;
    IdInterface2(any(Xi(gridHyPHMBasis.V0T(:, :)) == 1, 2) & any(Xi(gridHyPHMBasis.V0T(:, :)) == 3, 2)) = 1;

    temp = scatteredInterpolant(gridHyPHMBasis.baryT(:, 1), gridHyPHMBasis.baryT(:, 2), IdInterface1, 'linear');
    IdInterface1 = ceil(temp(gridHyPHM.baryT(:, 1), gridHyPHM.baryT(:, 2))) * k1;
    temp = scatteredInterpolant(gridHyPHMBasis.baryT(:, 1), gridHyPHMBasis.baryT(:, 2), IdInterface2, 'linear');
    IdInterface2 = ceil(temp(gridHyPHM.baryT(:, 1), gridHyPHM.baryT(:, 2))) * k2;


    for kT = find(IdInterface1 | IdInterface2)'
        L = 0;
        for i = 1:3
            if levelSetNew(gridHyPHM.V0E(gridHyPHM.E0T(kT, i), 1)) >= -eps & levelSetNew(gridHyPHM.V0E(gridHyPHM.E0T(kT, i), 2)) >= -eps
                L = L + gridHyPHM.areaE(gridHyPHM.E0T(kT, i)) / gridHyPHM.areaT(kT);
            end
        end
        SorceScale(kT) = L;
    end

    % nonlinear reaction equations

    epsilon2 = 10^(-18);
        %     nonlinearFunc = cell( 3, 1 );
    %     nonlinearFunc{1} = @(a,b,c) (IdInterface1.*(K1*max(b,sqrt(epsilon2))./max(a,sqrt(epsilon2))-1) ) .* SorceScale;
    %     nonlinearFunc{2} = @(a,b,c) (-IdInterface2.*( K2*b.*c-1) -IdInterface1.*(K1*b./max(a,sqrt(epsilon2))-1)).* SorceScale;
    %     nonlinearFunc{3} = @(a,b,c) (-IdInterface2.*(K2*b .* c-1) ) .* SorceScale;
    %
    %     nonlinearJacFunc = cell(3 , 3 );
    %     nonlinearJacFunc{1,1} = @(a,b,c) ((-IdInterface1.*(K1*b./max(a.^2,epsilon2))).*(a>sqrt(epsilon2)) + abs(sqrt(epsilon2)-a).* (a<sqrt(epsilon2))) .* SorceScale;
    %     nonlinearJacFunc{1,2} = @(a,b,c) IdInterface1.*(K1*1./max(a,sqrt(epsilon2))).* SorceScale;
    %     nonlinearJacFunc{1,3} = @(a,b,c) 0;
    %
    %     nonlinearJacFunc{2,1} = @(a,b,c) (IdInterface1.*(K1*b./max(a.^2,epsilon2))  ) .* SorceScale;
    %     nonlinearJacFunc{2,2} = @(a,b,c) ((-IdInterface2.*c *K2 ...
    %                                         -IdInterface1.*(K1./max(a,sqrt(epsilon2)))).* (b>sqrt(epsilon2))  + abs(sqrt(epsilon2)-b).* (b<sqrt(epsilon2))) .* SorceScale;
    %     nonlinearJacFunc{2,3} = @(a,b,c) -IdInterface2.* b .* K2 .* SorceScale;
    %
    %     nonlinearJacFunc{3,1} = @(a,b,c) 0;
    %     nonlinearJacFunc{3,2} = @(a,b,c) -IdInterface2 .* c .* K2  .* SorceScale;
    %     nonlinearJacFunc{3,3} = @(a,b,c) (-IdInterface2 .* b.* K2 ) .* SorceScale;

    nonlinearFunc = cell(3, 1);
    nonlinearFunc{1} = @(a, b, c) (IdInterface1 .* (K1 * max(b, sqrt(epsilon2)) ./ max(a, sqrt(epsilon2)) - 1)) .* SorceScale;
    nonlinearFunc{2} = @(a, b, c) (-IdInterface2 .* (K2 * b .* c - 1) - IdInterface1 .* (K1 * b ./ max(a, sqrt(epsilon2)) - 1)) .* SorceScale;
    nonlinearFunc{3} = @(a, b, c) (-IdInterface2 .* (K2 * b .* c - 1)) .* SorceScale;

    nonlinearJacFunc = cell(3, 3);
    nonlinearJacFunc{1, 1} = @(a, b, c) (-IdInterface1 .* (K1 * b ./ max(a.^2, epsilon2))) .* SorceScale;
    nonlinearJacFunc{1, 2} = @(a, b, c) IdInterface1 .* (K1 * 1 ./ max(a, sqrt(epsilon2))) .* SorceScale;
    nonlinearJacFunc{1, 3} = @(a, b, c) 0;

    nonlinearJacFunc{2, 1} = @(a, b, c) (IdInterface1 .* (K1 * b ./ max(a.^2, epsilon2))) .* SorceScale;
    nonlinearJacFunc{2, 2} = @(a, b, c) (-IdInterface2 .* c * K2 ...
        -IdInterface1 .* (K1 ./ max(a, sqrt(epsilon2)))) .* SorceScale;
    nonlinearJacFunc{2, 3} = @(a, b, c) -IdInterface2 .* b .* K2 .* SorceScale;

    nonlinearJacFunc{3, 1} = @(a, b, c) 0;
    nonlinearJacFunc{3, 2} = @(a, b, c) -IdInterface2 .* c .* K2 .* SorceScale;
    nonlinearJacFunc{3, 3} = @(a, b, c) (-IdInterface2 .* b .* K2) .* SorceScale;
    speciesCells = {ATransport; BTransport; CTransport};
    fprintf(['\n', 'CurrentTime: ', num2str(transportStepperSave.curtime), ' at step: ', num2str(transportStepperSave.curstep)])
    newtonIteration(speciesCells, nonlinearFunc, nonlinearJacFunc, 20);

    %save data
    AconcentrationSave.setdataGrid(transportStepperSave.curstep, ATransport.U.getdata(1), gridHyPHM);
    AFlux{transportStepperSave.curstep} = ATransport.Q.getdata(1);
    BconcentrationSave.setdataGrid(transportStepperSave.curstep, BTransport.U.getdata(1), gridHyPHM);
    BFlux{transportStepperSave.curstep} = BTransport.Q.getdata(1);
    CconcentrationSave.setdataGrid(transportStepperSave.curstep, CTransport.U.getdata(1), gridHyPHM);
    CFlux{transportStepperSave.curstep} = CTransport.Q.getdata(1);

    transportStepper.prev;

    if (currentTime >= finishTime - eps) | transportStepperSave.curstep == numtimeSteps
        break
    end

end

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

[interfaceLength, coordTripel] = evaluateInterface(cellGrid, Xi, distFunctions, false)

% add a circle of radius centering (X,Y) to the function handle for geometry evolution and
% adjust indicator Xi accordingly. Via two additional arguments the
% procedure generalizes to ellipses with ratio of half-axis and rotated by
% alpha
function [functionhandle, Xi] = drawcircle(X, Y, radius, epsilon, functionhandle, Xi, coord, varargin)
ratio = ones(size(X));
alpha = zeros(size(X));
if numel(varargin) == 1
    ratio = varargin{1};
    alpha = 0;
end
if numel(varargin) == 2
    ratio = varargin{1};
    alpha = varargin{2};
end
for i = 1:size(X, 2)
    rot = [cos(alpha(i)), sin(alpha(i)); -sin(alpha(i)), cos(alpha(i))];
    neu = @(x) min(radius(i)+epsilon-norm([ratio(i); 1] .* (rot * [x(1) - X(i); x(2) - Y(i)]), 2), -radius(i)+epsilon+norm([ratio(i); 1] .* (rot * [x(1) - X(i); x(2) - Y(i)]), 2));
    functionhandle = @(x) max(functionhandle(x), neu(x));
    Xi(radius(i)-sqrt((ratio(i)*rot(1, :)*(coord' - [X(i); Y(i)])).^2 + (rot(2, :) * (coord' - [X(i); Y(i)])).^2) > 0) = 3;
end
end
