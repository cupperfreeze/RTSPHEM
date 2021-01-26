function signedDistance = reinitializeLevelSet(grid, levelSet, isSigned, restrictDist, varargin)
%REINITIALIZELEVELSET    Reinitialize level set function.
%   Reinitializes a level set function LSF to the signed distance function
%   corresponding to the level set \{ LSF=0 \}. This SDF is the solution of the
%   Eikonal equation | \nabla U | = 1 with boundary conditions U = 0 on the
%   level set \{ LSF=0 \}.
%
%   This method uses the Fast Marching Method in combination with the
%   first-order Rouy-Tourin scheme for discretization of the gradient.
%
%   ATTENTION: This method only works on cartesian grids. Currently only
%   operational for positive Distances

assert(isa(grid, 'CartesianGrid'), ...
    'Grid is not of class ''CartesianGrid''.');
assert(abs(numel(levelSet) - grid.nodes) < eps, ...
    'Values are not compatible with the given grid.');

h = grid.stepSize;

% Assert that grid size is uniform in all dimensions.
assert(all(abs(h - mean(h)) < eps));

% Initialize values on/near the interface.
if length(varargin) == 0
    signedDistance = initializeFMM(grid, levelSet, isSigned);
else
    signedDistance = varargin{1};
end


% Initialize result vectors for the domains in which the level set function
% is positive / negative.
% Because the FMM is only defined for positive values of the solution
% (outward propagation), the array negativeDist saves the absolute value of
% the values in the negative domain of LSF. Then after the FMM is run, the
% minus sign is restored.
positiveDist = signedDistance;
negativeDist = -signedDistance;
isPositiveNode = (positiveDist > 0);
isNegativeNode = (negativeDist > 0);
positiveDist(isNegativeNode) = Inf;
negativeDist(isPositiveNode) = Inf;

% Mark all nodes directly on the interface as KNOWN.
isKnown = (abs(signedDistance) < eps(0));
grid.reshape(isKnown);
grid.reshape(signedDistance);

% Mark all nodes that are pre-initialized as TRIAL.
isTrialPos = false(numel(signedDistance), 1);
isTrialNeg = false(numel(signedDistance), 1);
isTrialPos(~isinf(positiveDist) & positiveDist > 0) = true;
isTrialNeg(~isinf(negativeDist) & negativeDist > 0) = true;


neighborIndices = [grid.pull{1, 1, 1}, grid.pull{2, 1, 1}, grid.pull{1, 1, 2}, grid.pull{2, 1, 2}]';

h = grid.stepSize;
dim = grid.dimension;

% Find the distances of the nodes on the positive/negative domain of LSF.
% use mex code if available
try
    [positiveDist] = ...
        FMM_iteration_mex(double(positiveDist), isKnown, ...
        h, uint32(dim), isPositiveNode, restrictDist, uint32(neighborIndices));
    [negativeDist] = ...
        FMM_iteration_mex(double(negativeDist), isKnown, ...
        h, uint32(dim), isNegativeNode, restrictDist, uint32(neighborIndices));
catch
    [positiveDist] = ...
        FMM_iteration(positiveDist, isKnown, ...
        h, uint32(dim), isPositiveNode, restrictDist, uint32(neighborIndices));
    [negativeDist] = ...
        FMM_iteration(negativeDist, isKnown, ...
        h, dim, isNegativeNode, restrictDist, uint32(neighborIndices));
end


signedDistance(isPositiveNode) = positiveDist(isPositiveNode);
% The sign of the right hand side determines whether the solution becomes
% the standard distance function or the signed distance function.
% For the standard distance: assign negativeDist( isNegativeNode ).
% For the signed distance: assign -negativeDist( isNegativeNode ).
signedDistance(isNegativeNode) = -negativeDist(isNegativeNode);

if (isa(grid, 'FoldedCartesianGrid'))
    signedDistance = grid.synchronizeValues(signedDistance);
end

end
