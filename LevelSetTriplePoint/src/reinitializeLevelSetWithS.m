function [Distance, orthogonalS] = reinitializeLevelSetWithS(grid, Distance, restrictDist, orthogonalS)
%REINITIALIZELEVELSET    Reinitialize level set function.
%   Reinitializes a level set function LSF to the signed distance function
%   corresponding to the level set \{ LSF=0 \}. This SDF is the solution of the
%   Eikonal equation | \nabla U | = 1 with boundary conditions U = 0 on the
%   level set \{ LSF=0 \}.
%
%   This method uses the Fast Marching Method in combination with the
%   first-order Rouy-Tourin scheme for discretization of the gradient.
%
%   ATTENTION: This method only works on cartesian grids.
%   orthogonalS is required to contain values at the nodes near the interface
%   levelSet is assumed to have noninf values at Interface nodes, elsewhere
%   inf

assert(isa(grid, 'CartesianGrid'), ...
    'Grid is not of class ''CartesianGrid''.');
assert(abs(numel(Distance) - grid.nodes) < eps, ...
    'Values are not compatible with the given grid.');

h = grid.stepSize;

% Assert that grid size is uniform in all dimensions.
assert(all(abs(h - mean(h)) < eps));


orthogonalS(isinf(Distance)) = inf;

% Initialize result vectors for the domains in which the level set function
% is positive / negative.
% Because the FMM is only defined for positive values of the solution
% (outward propagation), the array negativeDist saves the absolute value of
% the values in the negative domain of LSF. Then after the FMM is run, the
% minus sign is restored.
positiveDist = Distance;
isPositiveNode = (positiveDist > 0);

%
% Mark all nodes directly on the interface as KNOWN.
isKnown = (abs(Distance) < eps(0));
grid.reshape(isKnown);
grid.reshape(Distance);

% Mark all nodes that are pre-initialized as TRIAL.
isTrialPos = false(numel(Distance), 1);
isTrialPos(~isinf(positiveDist) & positiveDist > 0) = true;


neighborIndices = [grid.pull{1, 1, 1}, grid.pull{2, 1, 1}, grid.pull{1, 1, 2}, grid.pull{2, 1, 2}]';

h = grid.stepSize;
dim = grid.dimension;

%   Find the distances of the nodes on the positive domain of LSF.
%   use mex code if possible
try
    [positiveDist, orthogonalS] = ...
        FMM_iteration_WithS_mex(double(positiveDist), isKnown, ...
        h, uint32(dim), restrictDist, orthogonalS, uint32(neighborIndices));
catch
    [positiveDist, orthogonalS] = ...
        FMM_iteration_WithS(positiveDist, isKnown, ...
        h, uint32(dim), restrictDist, orthogonalS, uint32(neighborIndices));

end


Distance(isPositiveNode) = positiveDist(isPositiveNode);


if (isa(grid, 'FoldedCartesianGrid'))
    Distance = grid.synchronizeValues(Distance);
end

end
