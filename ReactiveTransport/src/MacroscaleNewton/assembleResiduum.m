function [residuum] = ...
    assembleResiduum(transportVariablesHyPHM, ...
    nonlinearFunc, Qpasts, Upasts, Qold, Uold, run)
%ASSEMBLERESIDUUM Assemble residuum for Newton solver
%   Generalizes the equations represented by class Transport of HyPHM to include
%   nonlinear source terms on the right hand side, i.e. equations of the form
%
%       \partial_t ( A u ) + \nabla \cdot q - F( u ) = R( q, u ) = 0,
%       q = - D \nabla u + C u + E
%
%   with initial resp. boundary conditions as in HyPHM. As F in general depends
%   on other variables that solve similar equations, this leads to a coupled
%   system of nonlinear equations which are solved via a global Newton method.
%
%   transportVariablesHyPHM: Cell array of variables of type Variable (from
%       numerical framework HyPHM).
%   nonlinearFunc: Cell array of function handles denoting the nonlinear
%       functions in F.
%
%   Returns:
%   - residuum: R evaluated with variables of the old time step.

numTransportVariables = numel(transportVariablesHyPHM);
unknownsCells = cell(numTransportVariables, 1);
nonlinearResiduumCells = cell(numTransportVariables, 1);
residuumCells = cell(numTransportVariables, 1);

[stiffnessDiagonalCells, loadVectorCells, degreesOfFreedomCells] = ...
    assembleHyPHMData(transportVariablesHyPHM, Upasts, run);

for i = 1:numTransportVariables

    transport = transportVariablesHyPHM{i};
    unknownsCells{i} = [Qold{i}; Uold{i}];
    nonlinearResiduumCells{i} = [zeros(transport.grid.numE, 1); ...
        transport.stepper.curtau * transport.grid.areaT .* nonlinearFunc{i}(Uold{:})];

    %         degreesOfFreedomCells{i} = ...
    %             logical( [ markFreeE; ones( transport.grid.numT, 1 ) ] );

    %         norm( stiffnessDiagonalCells{i}( degreesOfFreedomCells{i}, ...
    %                 degreesOfFreedomCells{i} ) * unknownsCells{i}( degreesOfFreedomCells{i} ) )
    %         norm( loadVectorCells{i}( degreesOfFreedomCells{i} ) )
    %         norm( nonlinearResiduumCells{i}( degreesOfFreedomCells{i} ) )

    residuumCells{i} = zeros(length(loadVectorCells{i}), 1);
    residuumCells{i}(degreesOfFreedomCells{i}) = ...
        stiffnessDiagonalCells{i}(degreesOfFreedomCells{i}, ...
        degreesOfFreedomCells{i}) ...
        * unknownsCells{i}(degreesOfFreedomCells{i}) ...
        -loadVectorCells{i}(degreesOfFreedomCells{i}) ...
        -nonlinearResiduumCells{i}(degreesOfFreedomCells{i});

end

residuum = cell2mat(residuumCells);

end
