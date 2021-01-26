%> @file +P0E/setfh2data.m Set data from f = @@(t, x), where @f$t@f$ [scalar], x @f$[2\times 1]@f$.
%> The dof's are on the barycenters of the triangles
function ret = setfh2data(time, grid, fun2)

ret = zeros(grid.numE, 1);
for kE = 1:grid.numE
    ret(kE) = fun2(time, grid.baryE(kE, :)');
end

end
