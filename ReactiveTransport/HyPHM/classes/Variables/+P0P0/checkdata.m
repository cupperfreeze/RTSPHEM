%> @file +P0P0/checkdata.m Check dimensions of raw numerical data.

function checkdata(grid, data)

msg = sprintf('HyPHM: Data has to have the dimension numT x 2, numT == %d.', grid.numT);
assert(isequal(size(data), [grid.numT, 2]), msg)

end
