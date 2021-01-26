%> @file +P2P2/checkdata.m Check dimensions of raw numerical data.
function checkdata(grid, data)

msg = sprintf(['HyPHM: Data has to have the dimension numV+numE x 2, ', ...
    'numV == %d, numE == %d.'], grid.numV, grid.numE);
assert(isequal(size(data), [grid.numV + grid.numE, 2]), msg)

end
