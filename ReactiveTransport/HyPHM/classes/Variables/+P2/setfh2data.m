%> @file +P2/setfh2data.m Set data from f = @@(t, x), where @f$t@f$ [scalar], x @f$[2\times 1]@f$.
function ret = setfh2data(time, grid, fun2)

ret = zeros(grid.numV+grid.numE, 1);
for kV = 1:grid.numV
    ret(kV) = fun2(time, grid.coordV(kV, :)');
end
for kE = 1:grid.numE
    ret(kE+grid.numV) = fun2(time, grid.baryE(kE, :)');
end

end
