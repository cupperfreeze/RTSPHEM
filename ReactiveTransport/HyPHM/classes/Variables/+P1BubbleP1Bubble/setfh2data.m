%> @file +P1BubbleP1Bubble/setfh2data.m Set data from f = @@(t, x), where @f$t@f$ [1\times 2], x @f$[2\times 1]@f$.
function ret = setfh2data(time, grid, fun2)


% DoFs are in the vertices and barycenters
valsOnV = zeros(grid.numV, 2);
for kV = 1:grid.numV
    valsOnV(kV, :) = fun2(time, grid.coordV(kV, :)');
end

valsOnT = zeros(grid.numT, 2);
for kT = 1:grid.numT
    valsOnT(kT, :) = fun2(time, grid.baryT(kT, :)') - mean(valsOnV(grid.V0T(kT, :)'));
end
ret = [valsOnV; valsOnT];

end
