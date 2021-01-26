%> @file +CR/checkdata.m Check dimensions of raw numerical data.
function checkdata(grid, data)

assert(isequal(size(data), [grid.numE, 1]), ...
    'HyPHM: Data has to have the dimension of numE x 1, numE == %d.', grid.numE)

end
