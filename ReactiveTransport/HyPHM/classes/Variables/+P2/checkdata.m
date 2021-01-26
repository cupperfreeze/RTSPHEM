%> @file +P2/checkdata.m Check dimensions of raw numerical data.
function checkdata(grid, data)

assert(isequal(size(data), [grid.numV + grid.numE, 1]), ...
    'HyPHM: Data has to have the dimension [numV+numE x 1], numV == %d, numE == %d.', ...
    grid.numV, grid.numE)

end
