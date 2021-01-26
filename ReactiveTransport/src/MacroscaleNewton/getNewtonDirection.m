function solutionVector = getNewtonDirection(transportVariablesHyPHM, ...
    jacobian, residuum, Qpasts, Upasts, Qold, Uold, run, iteration)
%SOLVENEWTONSTEP obtain direction for current Newton step

[~, ~, degreesOfFreedomCells, boundaryConditionCells, neumannCells] = ...
    assembleHyPHMData(transportVariablesHyPHM, Upasts, run);

degreesOfFreedom = cell2mat(degreesOfFreedomCells);

numRows = 0;
numTransportVariables = numel(transportVariablesHyPHM);
for i = 1:numTransportVariables

    numRows = numRows + size(boundaryConditionCells{i}, 1) ...
        +transportVariablesHyPHM{i}.grid.numT;

end
%     condest( jacobian( degreesOfFreedom, degreesOfFreedom ) )
solutionVector = zeros(numRows, 1);
%solutionVector2 = zeros( numRows, 1 );

%tic
%    solutionVector( degreesOfFreedom ) = ...
%        jacobian( degreesOfFreedom, degreesOfFreedom ) \ residuum( degreesOfFreedom );
%  toc


solutionVector(degreesOfFreedom) = solveNewton(jacobian(degreesOfFreedom, degreesOfFreedom), residuum(degreesOfFreedom), iteration, numTransportVariables);

%A=jacobian( degreesOfFreedom, degreesOfFreedom ) ;
%B=residuum( degreesOfFreedom );
%  norm(solutionVector-solutionVector2)
%  norm(solutionVector)
%    spy( jacobian( degreesOfFreedom, degreesOfFreedom ))
% return solutionVector

%     Qnew = cell( numTransportVariables, 1 );
%     Unew = cell( numTransportVariables, 1 );
%     startIndex = 1;
%     mask = false( length( solutionVector ), 1 );
%
%     Qdirection = cell( numTransportVariables, 1 );
%     Udirection = cell( numTransportVariables, 1 );
%
%     for i = 1:numTransportVariables
%
%         transport = transportVariablesHyPHM{i};
%         grid = transport.grid;
%         endIndexQ = startIndex + grid.numE - 1;
%
%         Qnew{i} = boundaryConditionCells{i};
%         mask( startIndex : endIndexQ ) = true;
%         dofQ = degreesOfFreedom & mask;
%         mask(:) = false;
%         dofQForCell = dofQ( startIndex : endIndexQ );
%         Qnew{i}( dofQForCell ) = Qold{i}( dofQForCell ) - solutionVector( dofQ );
%         Qdirection{i} = solutionVector( dofQ );
%
%         Unew{i} = zeros( grid.numT, 1 );
%         mask( ( endIndexQ + 1 ) : ( endIndexQ + grid.numT ) ) = true;
%         dofU = degreesOfFreedom & mask;
%         mask(:) = false;
%         dofUForCell = dofU( ( endIndexQ + 1 ) : ( endIndexQ + grid.numT ) );
%         Unew{i}( dofUForCell ) = Uold{i}( dofUForCell ) - solutionVector( dofU );
%         Udirection{i} = solutionVector( dofU );
%
%         startIndex = endIndexQ + grid.numT + 1;
%
%         st = transport.stepper;
%         dataC = transport.C.getdata( st.curstep );
%         dataE = transport.E.getdata( st.curstep );
%         Qnew{i}( neumannCells{i} ) = Qnew{i}( neumannCells{i} ) ...
%             + Unew{i}( grid.T0E( neumannCells{i}, 1 ) ) .* dataC( neumannCells{i} ) ...
%             + dataE( neumannCells{i} );
%
%     end


end
