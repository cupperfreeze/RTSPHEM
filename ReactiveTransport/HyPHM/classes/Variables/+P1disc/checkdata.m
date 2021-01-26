%> @file +P1disc/checkdata.m Check dimensions of raw numerical data, see Variable.type.
function checkdata(grid, data)

assert(isequal(size(data), [3 * grid.numT, 1]), ...
    'HyPHM: Data has to have the dimension (3*numT) x 1, (3*numT) == %d.', 3*grid.numT)

end
