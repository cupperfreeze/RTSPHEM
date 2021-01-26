function [Qnew, Unew] = applyNewtonStep(transportVariablesHyPHM, ...
    solutionVector, residuum, nonlinearFunc, Qpasts, Upasts, Qold, Uold, run)
%SOLVENEWTONSTEP Perform a single step of Newton's algorithm

[~, ~, degreesOfFreedomCells, boundaryConditionCells, neumannCells] = ...
    assembleHyPHMData(transportVariablesHyPHM, Upasts, run);

degreesOfFreedom = cell2mat(degreesOfFreedomCells);

numTransportVariables = numel(transportVariablesHyPHM);

%     startIndex = 1;
mask = false(length(solutionVector), 1);
armijoFactor = 1.0;
Qnew = cell(numTransportVariables, 1);
Unew = cell(numTransportVariables, 1);

while (true)
    startIndex = 1;

    for i = 1:numTransportVariables

        transport = transportVariablesHyPHM{i};
        grid = transport.grid;
        endIndexQ = startIndex + grid.numE - 1;

        Qnew{i} = boundaryConditionCells{i};
        mask(startIndex:endIndexQ) = true;
        dofQ = degreesOfFreedom & mask;
        mask(:) = false;
        dofQForCell = dofQ(startIndex:endIndexQ);
        Qnew{i}(dofQForCell) = Qold{i}(dofQForCell) ...
            -armijoFactor * solutionVector(dofQ);

        Unew{i} = zeros(grid.numT, 1);
        mask((endIndexQ+1):(endIndexQ + grid.numT)) = true;
        dofU = degreesOfFreedom & mask;
        mask(:) = false;
        dofUForCell = dofU((endIndexQ+1):(endIndexQ + grid.numT));
        Unew{i}(dofUForCell) = Uold{i}(dofUForCell) ...
            -armijoFactor * solutionVector(dofU);

        startIndex = endIndexQ + grid.numT + 1;

        st = transport.stepper;
        dataC = transport.C.getdata(st.curstep);
        dataE = transport.E.getdata(st.curstep);
        Qnew{i}(neumannCells{i}) = Qnew{i}(neumannCells{i}) ...
            +Unew{i}(grid.T0E(neumannCells{i}, 1)) .* dataC(neumannCells{i}) ...
            +dataE(neumannCells{i});

    end
    newResiduum = assembleResiduum(transportVariablesHyPHM, nonlinearFunc, ...
        Qpasts, Upasts, Qnew, Unew, run);

    if (norm(newResiduum) <= (1.0 - 1e-2 * armijoFactor) * norm(residuum))
        break;
    end
    armijoFactor = armijoFactor / 2;
    disp('        Armijo step size reduction.');

end

end
