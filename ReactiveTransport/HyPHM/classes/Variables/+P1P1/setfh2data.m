%> @file +P1P1/setfh2data.m Set data from f = @@(t, x), where @f$t@f$ [scalar], x @f$[2\times 1]@f$.
%> The dof's are on the verices and barycenters of the edges.
function ret = setfh2data(time, grid, fun2)

ret = zeros(grid.numV, 2);
for kV = 1:grid.numV
    ret(kV, 1:2) = fun2(time, grid.coordV(kV, :)');
end


end
