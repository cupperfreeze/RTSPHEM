%> @file +C/checkdata.m Check dimensions of raw numerical data.
function checkdata(~, data)

assert(isscalar(data), ...
    'HyPHM: Data has to be a scalar, but has dimensions %d x %d.', ...
    size(data, 1), size(data, 2))

end
