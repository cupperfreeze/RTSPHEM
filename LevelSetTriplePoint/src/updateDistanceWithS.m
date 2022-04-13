% similar to 'updateDistance', but additionally calculate velocity
% extension update
function [newDistance, NeworthogonalS] = updateDistanceWithS(h, dimensions, distance, ...
    updateIndex, orthogonalS, neighborIndices)

%UPDATEDISTANCE    Update distance values


dim = uint32(numel(dimensions));
dimd = double(dim);

minNeighbourVal = Inf(dim, numel(updateIndex));
DiscIndizes = Inf(dim, numel(updateIndex)); %Index in gradient discretization, 1=backward, 2=forward
backForthIndizes = cell(dim, 1); %Index of neighbor node
for i = 1:dim
    backForthIndizes{i} = Inf(2, numel(updateIndex));
    backForthIndizes{i}(1, :) = neighborIndices(2+2*(dimensions(i) - 1), updateIndex);
    backForthIndizes{i}(2, :) = neighborIndices(1+2*(dimensions(i) - 1), updateIndex);
end

for d = 1:dim

    for n = 1:numel(updateIndex)

        backwardIndex = neighborIndices(2+2*(dimensions(d) - 1), updateIndex(n));
        forwardIndex = neighborIndices(1+2*(dimensions(d) - 1), updateIndex(n));


        %backForthIndizes{d}(1,n) = backwardIndex;
        %backForthIndizes{d}(2,n) = forwardIndex;

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

        %calculate minimal neighbor distances and save where they
        %occured

        %[minNeighbourVal( d, n ), DiscIndizes( d, n )] = min( [distBackward, distForward] );
        if distBackward < distForward
            minNeighbourVal(d, n) = distBackward;
            DiscIndizes(d, n) = 1;
        else
            minNeighbourVal(d, n) = distForward;
            DiscIndizes(d, n) = 2;
        end

        % if both values are inf, this point will be disregarded later;
        % this is just to keep things running
        if isinf(distBackward) & isinf(distForward)
            temp = find([backwardIndex, forwardIndex] > 0.5, 1);
            DiscIndizes(d, n) = temp(1);
        end

    end

end


newDistance = distance;
NeworthogonalS = orthogonalS(updateIndex);


if (dim < 2) %for lower dimensional update (coearsend discretization)

    newDist = (minNeighbourVal + h(1))';
    newDistance = min(distance(updateIndex), newDist);
    GotUpdated = find(distance(updateIndex) > newDist);
    DiscIndexChosenX = backForthIndizes{1}(size(backForthIndizes{d}, 1) .* (0:length(updateIndex) - 1) + (DiscIndizes(1, :)));

    %In this approximtaion, S is constant along the edge
    %In case a node from the other interface side was chosen
    %to calculate distance  --> correct sign

    for p = 1:length(GotUpdated)
        if ~isinf(orthogonalS(DiscIndexChosenX(GotUpdated(p)))) & (abs(orthogonalS(updateIndex(GotUpdated(p))) + orthogonalS(DiscIndexChosenX(GotUpdated(p)))) < eps)
            NeworthogonalS(GotUpdated(p)) = -orthogonalS(DiscIndexChosenX(GotUpdated(p)));
        else
            NeworthogonalS(GotUpdated(p)) = orthogonalS(DiscIndexChosenX(GotUpdated(p)));
        end
    end

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

%Discretization Inidzes chosen for distance
DiscIndexChosenX = backForthIndizes{1}(size(backForthIndizes{1}, 1) .* (0:length(updateIndex) - 1) + DiscIndizes(1, :));
DiscIndexChosenY = backForthIndizes{2}(size(backForthIndizes{2}, 1) .* (0:length(updateIndex) - 1) + DiscIndizes(2, :));

denominator = ((minNeighbourVal(1, :) - newDist') ... %Calculate orthogonal S
    +(minNeighbourVal(2, :) - newDist'));


%In case a node from the other interface side was chosen
%to calculate distance  --> correct sign
signCorrectorX = ones(1, length(DiscIndexChosenX));
signCorrectorY = ones(1, length(DiscIndexChosenX));

for p = 1:length(DiscIndexChosenX)
    if ~isinf(orthogonalS(DiscIndexChosenX(p))) & (abs(orthogonalS(updateIndex(p)) + orthogonalS(DiscIndexChosenX(p))) < eps)
        signCorrectorX(p) = -1;
    end
    if ~isinf(orthogonalS(DiscIndexChosenY(p))) & (abs(orthogonalS(updateIndex(p)) + orthogonalS(DiscIndexChosenY(p))) < eps)
        signCorrectorY(p) = -1;
    end
end
%

newOrthS = (orthogonalS(DiscIndexChosenX) .* signCorrectorX .* (minNeighbourVal(1, :) - newDist') ...
    +orthogonalS(DiscIndexChosenY) .* signCorrectorY .* (minNeighbourVal(2, :) - newDist')) ...
    ./ denominator;


for p = 1:length(newOrthS) %If problems occur, e.g. vanishing gradient, assume S constant
    if isinf(newOrthS(p)) | abs(denominator(p)) < 10^(-8)
        newOrthS(p) = (orthogonalS(DiscIndexChosenX(p)) + orthogonalS(DiscIndexChosenY(p))) / 2;
    end
end

if (any(needsLowerUpdate)) %Perform lower dimensional update in each dimension, choose best

    potentialDistance = Inf(sum(needsLowerUpdate), dim);
    potentialorthogonalS = Inf(sum(needsLowerUpdate), dim);
    for d = 1:dim

        newDims = dimensions;
        newDims(d) = [];

        [potentialDistance(:, d), potentialorthogonalS(:, d)] = updateDistanceWithS(h, newDims, ...
            distance, updateIndex(needsLowerUpdate), orthogonalS, neighborIndices);

    end

    [newDist(needsLowerUpdate), DiscIndizes2] = ...
        min(potentialDistance(:, :), [], 2);

    newOrthS(needsLowerUpdate) = potentialorthogonalS(size(potentialorthogonalS, 1).*(DiscIndizes2 - 1)+(1:sum(needsLowerUpdate))');

end


% If lower dimensional update did improve, incorporate
newDistance = min(distance(updateIndex), newDist);
NeworthogonalS = newOrthS;


end