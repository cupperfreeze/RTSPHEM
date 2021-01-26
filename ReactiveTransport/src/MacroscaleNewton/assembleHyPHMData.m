function [stiffnessCells, loadVectorCells, degreesOfFreedomCells, ...
    boundaryConditionCells, neumannCells] = ...
    assembleHyPHMData(transportVariablesHyPHM, Upasts, run)
%ASSEMBLEHYPHMDATA assembly data from HyPHM transport instances


persistent a
persistent b
persistent c
persistent d
persistent e

%Perform assembly of static parts only for the first Newton step of a call
%and save them in persistent variables for further steps

if run == 0

    numTransportVariables = numel(transportVariablesHyPHM);

    stiffnessCells = cell(numTransportVariables, 1);
    loadVectorCells = cell(numTransportVariables, 1);
    degreesOfFreedomCells = cell(numTransportVariables, 1);
    neumannCells = cell(numTransportVariables, 1);
    boundaryConditionCells = cell(numTransportVariables, 1);


    for i = 1:numTransportVariables

        transport = transportVariablesHyPHM{i};
        stepper = transport.stepper;
        curTau = stepper.curtau;

        dataA = transport.A.getdata(stepper.curstep);
        dataAold = transport.A.getdata(stepper.curstep-1);
        dataAvold = NaN;
        dataB = transport.B.getdata(stepper.curstep);
        dataBold = transport.B.getdata(stepper.curstep-1);
        dataC = transport.C.getdata(stepper.curstep);
        dataD = transport.D.getdata(stepper.curstep);
        numT = transport.grid.numT;
        if strcmp(transport.D.type, 'CCCC') % if D is in CCCC, not P0P0P0P0
            dataD_P0P0P0P0 = zeros(numT, 4);
            row = reshape(dataD, 1, 4);
            for kT = 1:numT
                dataD_P0P0P0P0(kT, :) = row;
            end
            dataD = dataD_P0P0P0P0;
        else
            assert(strcmp(transport.D.type, 'P0P0P0P0'), 'HyPHM: Coefficient D has to have type CCCC or P0P0P0P0 but is %s.', transport.D.type)
        end
        assert(isequal(size(dataD), [numT, 4]), 'HyPHM kernel: Some strange things happen here.')

        if isequal(dataD, transport.D.getdata(stepper.curstep - 1))
            isDstationary = true;
        else
            isDstationary = false;
        end

        dataE = transport.E.getdata(stepper.curstep);
        dataF = transport.F.getdata(stepper.curstep);
        datauD = transport.uD.getdata(stepper.curstep);
        datagF = transport.gF.getdata(stepper.curstep);

        if (isa(transport, 'TransportLEVEL'))
            [markDirE, neumannCells{i}, markFluxE, markFreeE, ~, idxNeumE, ~, ~ markDirT] = ...
                transport.tp_indexsets();
            markNotDirT = ~markDirT;
        else
            [markDirE, neumannCells{i}, markFluxE, markFreeE, ~, idxNeumE, ~, ~] = ...
                transport.tp_indexsets();
            markNotDirT = ones(transport.grid.numT, 1);
        end

        isSlt = true;

        % Call HyPHM assembly routine for a single transport operator
        [B, C, D, E, bQ, bU, boundaryConditionCells{i}] = ...
            transport.assembleSystem( ...
            mat2cell(Upasts{i}, length(Upasts{i}), 1), ...
            1, curTau, ...
            dataA, dataAold, dataAvold, dataB, dataBold, dataC, dataD, ...
            isDstationary, dataE, dataF, datauD, datagF, markDirE, ...
            neumannCells{i}, markFluxE, idxNeumE, isSlt);

        stiffnessCells{i} = [B, C + D'; curTau * D, E];
        loadVectorCells{i} = [bQ; bU];

        degreesOfFreedomCells{i} = ...
            logical([markFreeE; markNotDirT]);


    end

    a = stiffnessCells;
    b = loadVectorCells;
    c = degreesOfFreedomCells;
    d = boundaryConditionCells;
    e = neumannCells;

else
    stiffnessCells = a;
    loadVectorCells = b;
    degreesOfFreedomCells = c;
    boundaryConditionCells = d;
    neumannCells = e;
end

end
