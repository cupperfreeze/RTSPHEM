%> @file +CR/setfh2data.m Set data from f = @@(t, x), where @f$t@f$ [scalar], x @f$[2\times 1]@f$.
function ret = setfh2data(time, grid, fun2)

% test the return value of fun2
assert(isscalar(fun2(time, grid.coordV(1, :)')), ...
    'HyPHM: The function handle for a Variable of type CR has to return a scalar.')

ret = zeros(grid.numE, 1);
for kE = 1:grid.numE
    ret(kE) = fun2(time, grid.baryE(kE, :)');
end

end
