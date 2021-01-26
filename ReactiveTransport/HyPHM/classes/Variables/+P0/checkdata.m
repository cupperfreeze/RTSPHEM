%> @file +P0/checkdata.m Check dimensions of raw numerical data.
function checkdata(grid, data)

assert(isequal(size(data), [grid.numT, 1]), ...
    'HyPHM: Data has to have the dimension numT x 1, numT == %d.', grid.numT)

end
