%> @file +P0/setfh2data.m Set data from f = @@(t, x), where @f$t@f$ [scalar], x @f$[2\times 1]@f$.
%> The dof's are on the barycenters of the triangles
function ret = setfh2data(time, grid, fun2)

accelerate=false;
try
    accelerate = isequal(size(fun2(time, grid.baryT(1:2, :)')), [2, 1]);
end
if accelerate
    ret = fun2(time, grid.baryT'); %accelerate if vectorized
else
    ret = zeros(grid.numT, 1);
    for kT = 1:grid.numT
        ret(kT) = fun2(time, grid.baryT(kT, :)');
    end
end

end
