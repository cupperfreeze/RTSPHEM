function [transportVariablesHyPHM, iteration] = newtonIteration(transportVariablesHyPHM, ...
    nonlinearFunc, nonlinearJacobianFunc, numIterations)
%Solve system of 'numTransportVariables' transport equations coupled by a
%nonlinear source term. Routine expects cell input of transport instances

numTransportVariables = numel(transportVariablesHyPHM);

Qpasts = cell(numTransportVariables, 1);
Upasts = cell(numTransportVariables, 1);
Q = cell(numTransportVariables, 1);
U = cell(numTransportVariables, 1);
for i = 1:numTransportVariables
    transport = transportVariablesHyPHM{i};
    Qpasts{i} = transport.Q.getdata(transport.stepper.curstep-1);
    Upasts{i} = transport.U.getdata(transport.stepper.curstep-1);
    Q{i} = Qpasts{i};
    U{i} = Upasts{i};
end

iteration = 0;
isConverged = false;
relativeTolerance = 1e-6;
absoluteTolerance = 1e-13;
fprintf('\n');

run = 0; % control if matrix assembly is already peformed --> use stored values
defaultResiduum = assembleResiduum(transportVariablesHyPHM, nonlinearFunc, ...
    Qpasts, Upasts, Q, U, run);
run = 1;
while (iteration < numIterations && ~isConverged)

    residuum = assembleResiduum(transportVariablesHyPHM, nonlinearFunc, ...
        Qpasts, Upasts, Q, U, run);

    fprintf('Newton step %i: Residuum = %d\n', iteration, norm(residuum));
    if (norm(residuum) < relativeTolerance * norm(defaultResiduum) ...
            +absoluteTolerance)
        %             isConverged = true;
        break;
    end

    jacobian = assembleJacobian(transportVariablesHyPHM, ...
        nonlinearJacobianFunc, Qpasts, Upasts, Q, U, run);

    %         [ Qnew, Unew ] = getNewtonDirection( transportVariablesHyPHM, ...
    %             jacobian, residuum, Qpasts, Upasts, Q, U );
    solutionVector = getNewtonDirection( ...
        transportVariablesHyPHM, jacobian, residuum, Qpasts, Upasts, Q, U, run, iteration);

    [Qnew, Unew] = applyNewtonStep(transportVariablesHyPHM, ...
        solutionVector, residuum, nonlinearFunc, Qpasts, Upasts, Q, U, run);
    %         isConverged = true;
    %         for i = 1:numTransportVariables
    %
    %             fprintf( '    Distance of flux iterates: %d\n', norm( Qnew{i} - Qpasts{i} ) );
    %             if ( norm( Qnew{i} - Q{i} ) > tolerance )
    %                 isConverged = false;
    %             end
    %
    %             fprintf( '    Distance of scalar iterates: %d\n\n', norm( Unew{i} - Upasts{i} ) );
    %             if ( norm( Unew{i} - U{i} ) > tolerance )
    %                 isConverged = false;
    %                 break;
    %             end
    %
    %         end
    iteration = iteration + 1;

    for i = 1:numTransportVariables
        Q{i} = Qnew{i};
        U{i} = Unew{i};
    end
end

if iteration > 0 % clear last preconditioner
    solveNewton(sparse(1), 0, 0, 0, 0);
end
clear solveNewton

for i = 1:numTransportVariables

    transport = transportVariablesHyPHM{i};
    transport.Q.setdata(transport.stepper.curstep, Q{i});
    transport.U.setdata(transport.stepper.curstep, max(U{i}, eps)); % concentrations regularization

end

end
