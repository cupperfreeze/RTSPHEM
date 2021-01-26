%> @file +P1disc/setfh2data.m Set data from f = @@(t, x), where @f$t@f$ [scalar], x @f$[2\times 1]@f$.
function ret = setfh2data(time, grid, fun2)

% test the return value of fun2
assert(isscalar(fun2(time, grid.coordV(1, :)')), ...
    'HyPHM: The function handle for a Variable of type P1 has to return a scalar.')

% although elements of P1disc have DOFs on each vertex FOR EACH TRIANGLE
% that are not shared, the value of a fixed vertex is the very same for
% each surrounding triangle due to the fact that fun2 returns a single
% value.
valsOnV = zeros(grid.numV, 1);
for kV = 1:grid.numV
    valsOnV(kV) = fun2(time, grid.coordV(kV, :)');
end
ret = reshape(valsOnV(grid.V0T)', 3*grid.numT, 1);

end
