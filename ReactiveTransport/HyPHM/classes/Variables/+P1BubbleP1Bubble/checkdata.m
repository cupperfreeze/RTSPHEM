%> @file +P1BubbleP1Bubble/checkdata.m Check dimensions of raw numerical data, see Variable.type.
function checkdata(grid, data)

assert(isequal(size(data), [grid.numV + grid.numT, 2]), ...
    'HyPHM: Data has to have the dimension (numV+numT) x 2, (numV+numT) == %d.', grid.numV+grid.numT)

end
