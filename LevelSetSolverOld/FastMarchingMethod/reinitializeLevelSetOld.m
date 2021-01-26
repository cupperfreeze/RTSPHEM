function signedDistance = reinitializeLevelSet(grid, levelSet)
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

assert(isa(grid, 'CartesianGrid'), ...
    'Grid is not of class ''CartesianGrid''.');
assert(abs(numel(levelSet) - grid.nodes) < eps, ...
    'Values are not compatible with the given grid.');

% Initialize values on/near the interface.
signedDistance = initializeFMM(grid, levelSet);

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
trialIndexPos = find(isTrialPos);
trialIndexNeg = find(isTrialNeg);

% Find the distances of the nodes on the positive domain of LSF.
while (any(isTrialPos(:)))
    [positiveDist, isKnown, isTrialPos, trialIndexPos] = ...
        FMM_iteration(positiveDist, isKnown, isTrialPos, ...
        trialIndexPos, grid, isPositiveNode);

    %         % Find the TRIAL node with the minimal distance and set it to KNOWN.
    %         [~, minIndex] = min( positiveDist( isTrialPos ) );
    %         curIndex = trialIndexPos( minIndex );
    %         isKnown( curIndex ) = true;
    %         isTrialPos( curIndex ) = false;
    %
    %         % Mark all neighbours of the newly KNOWN node and mark them as TRIAL
    %         % (if they are not yet KNOWN).
    %         newIndices = NaN( 2, grid.getDimension );
    %         for d = 1:grid.getDimension
    %             newIndices(2*d-1) = grid.getBackward( curIndex, d );
    %             newIndices(2*d) = grid.getForward( curIndex, d );
    %         end
    %         isTrialPos( newIndices ) = ~isnan( newIndices ) ...
    %             & ~isKnown( newIndices ) ...
    %             & isPositiveNode( newIndices );
    %         trialIndexPos = find( isTrialPos );
    %
    %         % Update the distance values for all new TRIAL nodes.
    %         positiveDist = updateDistance( grid, 1:grid.getDimension, ...
    %             positiveDist, newIndices(:) );
end

% Find the distances of the points on the negative domain of LSF.
while (any(isTrialNeg(:)))
    [negativeDist, isKnown, isTrialNeg, trialIndexNeg] = ...
        FMM_iteration(negativeDist, isKnown, isTrialNeg, ...
        trialIndexNeg, grid, isNegativeNode);

    %         [~, maxIndex] = max( negativeDist( isTrialNeg ) );
    %         curIndex = trialIndexNeg( maxIndex );
    %         isKnown( curIndex ) = true;
    %         isTrialNeg( curIndex ) = false;
    %
    %         newIndices = NaN( 2, grid.getDimension );
    %         for d = 1:grid.getDimension
    %             newIndices(2*d-1) = grid.getBackward( curIndex, d );
    %             newIndices(2*d) = grid.getForward( curIndex, d );
    %         end
    %         isTrialNeg( newIndices ) = ~isnan( newIndices ) ...
    %             & ~isKnown( newIndices ) ...
    %             & isNegativeNode( newIndices );
    %         trialIndexNeg = find( isTrialNeg );
    %
    %         negativeDist = -updateDistance( grid, 1:grid.getDimension, ...
    %             -negativeDist, newIndices(:) );
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

function [distance, isKnown, isTrial, trialIndex] = ...
    FMM_iteration(distance, isKnown, isTrial, trialIndex, grid, toBeComputed)

% Find the TRIAL node with the minimal distance and set it to KNOWN.
[~, minIndex] = min(distance(isTrial));
curIndex = trialIndex(minIndex);
isKnown(curIndex) = true;
isTrial(curIndex) = false;

% Mark all neighbours of the newly KNOWN node and mark them as TRIAL
% (if they are not yet KNOWN).
newIndices = NaN(2, grid.dimension);
for d = 1:grid.dimension
    newIndices(2*d-1) = grid.getBackward(curIndex, d);
    newIndices(2*d) = grid.getForward(curIndex, d);
end
newIndices(isnan(newIndices)) = [];
isTrial(newIndices(~isnan(newIndices))) = ...
    ~isKnown(newIndices) ...
    & toBeComputed(newIndices);
trialIndex = find(isTrial);

% Update the distance values for all new TRIAL nodes.
distance = updateDistance(grid, 1:grid.dimension, ...
    distance, newIndices(:));

end
