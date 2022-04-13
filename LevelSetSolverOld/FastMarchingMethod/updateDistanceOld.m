function [newDistance] = updateDistanceOld(grid, dimensions, distance, ...
    updateIndex)
%UPDATEDISTANCE    Update distance values

dim = numel(dimensions);
h = grid.stepSize;

% Assert that grid size is uniform in all dimensions.
assert(all(abs(h - mean(h)) < eps));

minNeighbourVal = Inf(dim, numel(updateIndex));

for d = 1:dim

    for n = 1:numel(updateIndex)

        backwardIndex = grid.getBackward(updateIndex(n), dimensions(d));
        forwardIndex = grid.getForward(updateIndex(n), dimensions(d));

        if (isnan(backwardIndex))
            distBackward = Inf;
        else
            distBackward = distance(backwardIndex);
        end
        if (isnan(forwardIndex))
            distForward = Inf;
        else
            distForward = distance(forwardIndex);
        end

        minNeighbourVal(d, n) = min(distBackward, distForward);

    end

end


newDistance = distance;

if (dim < 2)

    newDist = (minNeighbourVal + h(1))';
    newDistance(updateIndex) = min(distance(updateIndex), newDist);

    return;

end

%     newDist = min( minNeighbourVal + repmat( grid.getH', 1, ...
%         numel(updateIndex) ) )';

l1Norms = sum(minNeighbourVal)';
l2Norms = sqrt(sum(minNeighbourVal.^2))';
lInfNorms = max(minNeighbourVal)';

dim = double(dim);
discriminant = dim * h(1)^2 - dim * l2Norms.^2 + l1Norms.^2;
isSolvable = (discriminant >= 0);
isConsistent = (h(1)^2 >= dim * lInfNorms.^2 + l2Norms.^2 ...
    -2 * l1Norms .* lInfNorms);
needsLowerUpdate = ~isSolvable | ~isConsistent;

newDist = (l1Norms + sqrt(discriminant)) / dim;

if (any(needsLowerUpdate))

    potentialDistance = Inf(numel(distance), dim);

    for d = 1:dim

        newDims = dimensions;
        newDims(d) = [];

        potentialDistance(:, d) = updateDistanceOld(grid, newDims, ...
            distance, updateIndex(needsLowerUpdate));

    end

    newDist(needsLowerUpdate) = ...
        min(potentialDistance(updateIndex(needsLowerUpdate), :), [], 2);

end

altDist = newDist;

altDist(needsLowerUpdate) = min(minNeighbourVal(:, needsLowerUpdate)+ ...
    repmat(h', 1, sum(needsLowerUpdate)))';

%     needsLowerUpdate
%     [newDist, altDist, newDist - altDist]

newDistance(updateIndex) = min(distance(updateIndex), newDist);

end
