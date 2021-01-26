%> @file +CCCC/checkdata.m Check dimensions of raw numerical data.

function checkdata(~, data)

assert(isequal(size(data), [2, 2]), ...
    'HyPHM: Data has to be a 2x2 matrix, but has dimensions %d x %d.', size(data, 1), size(data, 2))

end
