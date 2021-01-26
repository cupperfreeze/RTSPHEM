%> @file +P1/checkdata.m Check dimensions of raw numerical data.
function checkdata(grid, data)

assert(isequal(size(data), [grid.numV, 1]), ...
    'HyPHM: Data has to have the dimension numV x 1, numV == %d.', grid.numV)

end
