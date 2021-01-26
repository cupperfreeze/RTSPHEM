%> @file +P0P0P0P0/setfh2data.m Set data from f = @@(t, x), where @f$t@f$ [scalar], x @f$[2\times 1]@f$.
%> The dof's are on the barycenters of the triangles
function ret = setfh2data(time, grid, fun2)

% check return value of fh
msg = 'HyPHM: The function_handle used to set the data of a Variable of type P0P0P0P0 has to return a 2x2-matrix.';
%assert(isequal(size(fun2(time, grid.baryT(1, :)')), [2,2]), msg)

if ~isequal(size(fun2(time, grid.baryT(1, :)')), [2, 2]) %accelerate if vectorized
    ret = fun2(time, grid.baryT');
else
    ret = zeros(grid.numT, 4);
    for kT = 1:grid.numT
        ret(kT, 1:4) = reshape(fun2(time, grid.baryT(kT, :)'), 1, 4);
    end

end
end
