%> @file checkdir.m Checks if the call of HyPHM was made from the root directory.

w = what;
if ~any(ismember(w.m, 'HyPHM.m'))
    s1 = 'HyPHM: HyPHM has to be executed in its root directory.  Change to ';
    s2 = which('HyPHM.m');
    s2 = s2(1:end-length('/HyPHM.m'));
    error([s1, s2])
end
