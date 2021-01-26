%>  @file +P1P1/checkdata.m Check dimensions of raw numerical data.
function ret = checkdata(grid, data)

msg = sprintf(['HyPHM: Data has to have the dimension numV x 2, ', ...
    'numV == %d, however, size(data) = [%d, %d].'], grid.numV, size(data, 1), size(data, 2));
assert(isequal(size(data), [grid.numV, 2]), msg)

end
