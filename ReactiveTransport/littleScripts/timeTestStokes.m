%script to time test DolomiteMicroMacro and save resulting calcium concentration
%requires

%timeStorage = NaN(20,7);
global timeNewton;
global timeLevelSet;
global timeDiff;
global timePerm;
global timeDarcy;
global hydroData;
global numPartitionsMicroscale;
global macroscaleStepSize;
global PermeabilityTensorTest;
global CalciumConcentrationTest;
global GridTest;
global numProc;

numProc = 1;

% Different microscopic and macroscopic mesh parameter combintaions
Micro = kron([25, 50, 100, 200, 400], [1, 1, 1, 1]);
Macro = repmat([1 / 150, 1 / 300, 1 / 600, 1 / 1200], 1, 5);

for me = 1:numel(Micro)
    % numProc = 2^(me-1);
    numPartitionsMicroscale = Micro(me);
    macroscaleStepSize = Macro(me);
    timeNewton = 0;
    timeLevelSet = 0;
    timeDiff = 0;
    timePerm = 0;
    timeDarcy = 0;
    hydroData = 0;
    tot = tic;
    DolomiteMicroMacro;
    timeTotal = toc(tot);

    timeStorage(me, 1) = timeNewton;
    timeStorage(me, 2) = timeLevelSet;
    timeStorage(me, 3) = timeDiff;
    timeStorage(me, 4) = timePerm;
    timeStorage(me, 5) = timeDarcy;
    timeStorage(me, 6) = timeTotal;
    timeStorage(me, 7) = hydroData;
    PermeabilityTensors{me} = PermeabilityTensorTest;
    CalciumConcentration{me} = CalciumConcentrationTest;
    Grids{me} = GridTest;

end
