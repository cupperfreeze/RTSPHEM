%> @file +P1/setfh2data.m Set data from f = @@(t, x), where @f$t@f$ [scalar], x @f$[2\times 1]@f$.
function ret = setfh2data(time, grid, fun2)

assert(isscalar(fun2(time, grid.coordV(1, :)')), ...
    'HyPHM: The function handle for a Variable of type P1 has to return a scalar.')

ret = zeros(grid.numV, 1);
for kV = 1:grid.numV
    ret(kV) = fun2(time, grid.coordV(kV, :)');
end

end
