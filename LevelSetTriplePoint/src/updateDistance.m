function [newDistance] = updateDistance(h, dimensions, distance, ...
    updateIndex, neighborIndices)
%UPDATEDISTANCE    Update distance values in points specified by updateIndex

dim = uint32(numel(dimensions));
dimd = double(dim);

minNeighbourVal = Inf(dim, numel(updateIndex));

for d = 1:dim

    for n = 1:numel(updateIndex)

        backwardIndex = neighborIndices(2+2*(dimensions(d) - 1), updateIndex(n));
        forwardIndex = neighborIndices(1+2*(dimensions(d) - 1), updateIndex(n));

        if (backwardIndex < 0.5)
            distBackward = Inf;
        else
            distBackward = distance(backwardIndex);
        end
        if (forwardIndex < 0.5)
            distForward = Inf;
        else
            distForward = distance(forwardIndex);
        end

        minNeighbourVal(d, n) = min(distBackward, distForward);

    end

end


newDistance = distance(updateIndex);

% Switch to coarser discretization of gradient
if (dim < 2)

    newDist = (minNeighbourVal + h(1))';
    newDistance = min(distance(updateIndex), newDist);

    return;

end

%     newDist = min( minNeighbourVal + repmat( grid.getH', 1, ...
%         numel(updateIndex) ) )';

l1Norms = sum(minNeighbourVal)';
l2Norms = sqrt(sum(minNeighbourVal.^2))';
lInfNorms = max(minNeighbourVal)';


discriminant = dimd * h(1)^2 - dimd * l2Norms.^2 + l1Norms.^2;
isSolvable = (discriminant >= 0);
isConsistent = (h(1)^2 >= dimd * lInfNorms.^2 + l2Norms.^2 ...
    -2 * l1Norms .* lInfNorms);
needsLowerUpdate = ~isSolvable | ~isConsistent;

newDist = (l1Norms + sqrt(discriminant .* isSolvable)) / dimd;

% first order discretization did not succeed, coarsen discretization of gradient
if (any(needsLowerUpdate))

    potentialDistance = Inf(sum(needsLowerUpdate), dim);

    for d = 1:dim

        newDims = dimensions;
        newDims(d) = [];

        potentialDistance(:, d) = updateDistance(h, newDims, ...
            distance, updateIndex(needsLowerUpdate), neighborIndices);

    end

    newDist(needsLowerUpdate) = ...
        min(potentialDistance, [], 2);

end


newDistance = min(distance(updateIndex), newDist);

end
