% Procedure to evolve unit cell geometries in one macroscopic timestep

function microCell = evolveStep(microCell, interfaceVelocities, ...
    timeStep, curtau, numPartitionsMicroscale)

% adjust timestep = position in memory for current data
if microCell.reduceData
    timeStep = 1;
else
    microCell.memoryStep = microCell.memoryStep + 1;
end

restrictDist = microCell.restrictDist;
numberOfSlices = numel(microCell.Xi);

% create local storage
localreInit = cell(numberOfSlices, 1);
localPhi = cell(numberOfSlices, 1);
localXi = cell(numberOfSlices, 1);
localsigned = cell(numberOfSlices, 1);
localadapS = cell(numberOfSlices, 1);
localadapTime = cell(numberOfSlices, 1);
localinterfaceLength = cell(numberOfSlices, 1);
localdistFunctions = cell(numberOfSlices, 1);
CFL = inf(numberOfSlices, 1);
AdapCounter = zeros(numberOfSlices, 1); %mark which cells were evolved

parfor index = 1:numberOfSlices

    currentTime = 0;
    %unpack variables
    Phi = microCell.Phi{index}(:, timeStep);
    Xi = microCell.Xi{index}(:, timeStep);
    Lev = microCell.signed{index}(:, timeStep);
    cellGrid = microCell.cellGrid;
    distFunctions = microCell.distFunctions{index};
    reInit = microCell.reInit(index);
    interfaceLength = microCell.interfaceLength{index};

    %Adaptivity scheme
    localadapS{index} = interfaceVelocities{index} * curtau + microCell.adapS{index};
    localadapTime{index} = microCell.adapTime{index} + curtau;

    if max(localadapS{index}(:)) > 1 / numPartitionsMicroscale / 5

        AdapCounter(index) = 1;
        T = localadapTime{index};
        interfaceVelocities{index} = localadapS{index} / T;
        CFL(index) = 1 * 1 / 3 * 1 / numPartitionsMicroscale / max(interfaceVelocities{index}(:));
        localadapS{index} = zeros(3);
        localadapTime{index} = 0;

        % VIIM step
        for j = 1:ceil(T/CFL(index))

            StepSizeTime = min(T-currentTime, CFL(index));
            currentTime = currentTime + StepSizeTime;

            %Reconstruct Vonoroi Interface and velocity extension
            [dV, S] = step234(distFunctions, cellGrid, Xi, interfaceVelocities{index}, 8/cellGrid.nodesPerDimension(1));

            %Update Phi via dV, but only every 10 steps
            if mod(reInit, 10) == 0
                Phi = microCell.epsilon - dV;
            end
            reInit = reInit + 1;

            %Evolve Phi
            Phi = levelSetEquationTimeStep(StepSizeTime, 0, Phi, ...
                cellGrid, S', 1);


            %Update Xi according to Phi,
            %Calculate distance function from each subdomain
            [distFunctions, Xi] = step1(Xi, Phi, cellGrid, restrictDist);

        end

        [interfaceLength, ~] = evaluateInterface(cellGrid, Xi, distFunctions, true);

    end
  
    localinterfaceLength{index} = interfaceLength(1:2);
    localreInit{index} = reInit;
    localPhi{index} = Phi;
    localXi{index} = Xi;
    localdistFunctions{index} = distFunctions;

    if ~isinf(CFL(index))
        Lev = dV; %calculate approximation to signed distance function w.r.t solid
        Lev(Xi == 1) = -Lev(Xi == 1);
        Lev(Xi == 2) = max(Lev(Xi == 2), 10^(-4));
        Lev(Xi == 3) = max(Lev(Xi == 3), 10^(-4));
    end

    localsigned{index} = Lev;

end

if ~microCell.reduceData
    timeStep = timeStep + 1;
end

% pack
for index = 1:numberOfSlices
    microCell.adapS{index} = localadapS{index};
    microCell.adapTime{index} = localadapTime{index};
    microCell.reInit(index) = localreInit{index};
    microCell.Phi{index}(:, timeStep) = localPhi{index};
    microCell.Xi{index}(:, timeStep) = localXi{index};
    microCell.signed{index}(:, timeStep) = localsigned{index};
    microCell.interfaceLength{index}(timeStep, :) = localinterfaceLength{index};
    microCell.distFunctions{index} = single(localdistFunctions{index});
end
disp(['LevelSet Adaptivity: ', ...
    num2str(numberOfSlices - sum(AdapCounter)), ' out of ', num2str(numberOfSlices), ' FMM saved']);
microCell.saved = {microCell.saved{:}, AdapCounter};
end
