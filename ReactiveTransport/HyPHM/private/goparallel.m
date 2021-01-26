%> @file goparallel.m Initialization of Matlab-workers, if admissible.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> Initialization of Matlab-workers, if admissible.

printline(1, 'Initialization of Matlab-Workers');

try
    if matlabpool('size') == 0
        fprintf('     '); % text align space
        if exist('maxworkers', 'var')
            matlabpool('open', 'local', num2str(maxworkers)) % open new pool with desired number of local workers
        else
            matlabpool('open', 'local') % open new pool with max number of local workers
        end
        printline(3, 'Number of available cores: %d', matlabpool('size'))
    end
    if matlabpool('size') < 12
        printline(2, 'You may want to set the number of maximum workers (up to 12!) by:')
        printline(3, '>> myCluster = parcluster(''local'');')
        printline(3, '>> myCluster.NumWorkers = 12')
        printline(3, '>> saveProfile(myCluster);\n')
    end
catch
    printline(-1, 'HyPHM: No parallelization possible for this Matlab/Octave version!')
    printline(3, 'Check HyPHM/private/goparallel if you want to go into detail here.')
end
