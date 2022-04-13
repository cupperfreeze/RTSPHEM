function [jacobian] = ...
    assembleJacobian(transportVariablesHyPHM, ...
    nonlinearJacobianFunc, Qpasts, Upasts, Qold, Uold, run)
%   ASSEMBLEJACOBIAN Assemble jacobian and residuum for Newton solver
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
%   nonlinearJacobianFunc: Cell matrix of function handles denoting the partial
%       derivatives of the nonlinear functions in F. Function handles must take
%       #transportVariables arguments and must be able to take vector-sized
%       inputs. Size: #transportVariables x #transportVariables
%
%   Returns:
%   - jacobian: Jacobian of R evaluated with variables of the old time step.
%   - residuum: R evaluated with variables of the old time step.

numTransportVariables = numel(transportVariablesHyPHM);

[stiffnessCells, ~, ~] = ...
    assembleHyPHMData(transportVariablesHyPHM, Upasts, run);
jacobianNonlinearPart = cell(numTransportVariables, numTransportVariables);

for i = 1:numTransportVariables

    transport = transportVariablesHyPHM{i};

    %         [ B, C, D, E, bQ, bU, boundaryConditions{i} ] = ...
    %             transport.assembleSystem( ...
    %                 mat2cell( Uold{i}, length( Uold{i} ), 1 ), ...
    %                 1, curTau, ...
    %                 dataA, dataAold, dataAvold, dataB, dataBold, dataC, dataD, ...
    %                 isDstationary, dataE, dataF, datauD, datagF, markDirE, ...
    %                 markNeumE{i}, markFluxE, idxNeumE, isSlt );

    indices = (1:transport.grid.numT) + transport.grid.numE;

    % obtain linearized source term
    for j = 1:numTransportVariables
        assert(transport.grid.numT == transportVariablesHyPHM{j}.grid.numT);
        jacobianNonlinearPart{i, j} = sparse(indices, indices, ...
            transport.stepper.curtau*transport.grid.areaT.*nonlinearJacobianFunc{i, j}(Uold{:}));

    end

end

%     jacobian = kronVec( eye( numTransportVariables ), ...
%         cell2mat( stiffnessCells ) ) - cell2mat( jacobianNonlinearPart );
jacobian = blkdiag(stiffnessCells{:}) - cell2mat(jacobianNonlinearPart);

end
