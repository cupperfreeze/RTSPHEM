function [distance, orthogonalS] = ...
    FMM_iteration_WithS(distance, isKnown, h, dim, restrictDist, orthogonalS, neighborIndices)

%%   Initilaize heap

nonZeroDist = distance;
nonZeroDist(isKnown) = inf;
[val, ind] = sort(nonZeroDist);
numelements = numel(val);
sizeOfHeap = round(2^(ceil(log2(numelements + 1)))) - 1;

valuesHeap = inf(sizeOfHeap, 1);
valuesHeap(1:numelements) = val;

indexHeap = uint32(inf(sizeOfHeap, 1));
indexHeap(1:numelements) = uint32(ind);

whereIsWhat = uint32(inf(sizeOfHeap, 1));
whereIsWhat(indexHeap(1:numelements)) = uint32(1:numelements);

kHeap = uint32(sum(val < inf));

doExit = false;

%%

while (kHeap > 0)
    % Find the TRIAL node with the minimal distance and set it to KNOWN.
    %   [valueExit, minIndex] = min( distance( isTrial ) );
    %   find(isTrial);

    %% Extract Min

    GMMconst = 200;
    HEAPMin = zeros(1, GMMconst);
    HEAPMinIndex = zeros(1, GMMconst);

    for m = 1:GMMconst
        HEAPMin(m) = valuesHeap(1);
        HEAPMinIndex(m) = indexHeap(1);

        if HEAPMin(m) > HEAPMin(1) + 1 / sqrt(2) * h(1) / 1.01 || HEAPMin(m) > restrictDist
            if ~isinf(HEAPMin(m)) && HEAPMin(m) > restrictDist
                doExit = true;
            end
            HEAPMin = HEAPMin(1:(m - 1));
            HEAPMinIndex = HEAPMinIndex(1:(m - 1));
            break
        end


        valuesHeap(1) = inf;
        kHeap = kHeap - 1;
        i = uint32(1);

        while i < double(sizeOfHeap) / 2


            if valuesHeap(2*i) <= valuesHeap(2*i+1) && valuesHeap(i) > valuesHeap(2*i)
                whereIsWhat(indexHeap(2 * i)) = i;
                whereIsWhat(indexHeap(i)) = 2 * i;

                swap = valuesHeap(2*i);
                valuesHeap(2*i) = valuesHeap(i);
                valuesHeap(i) = swap;

                swap = indexHeap(2*i);
                indexHeap(2*i) = indexHeap(i);
                indexHeap(i) = swap;

                i = 2 * i;
                continue
            else

                if valuesHeap(2*i) >= valuesHeap(2*i+1) && valuesHeap(i) > valuesHeap(2*i+1)
                    whereIsWhat(indexHeap(2 * i + 1)) = i;
                    whereIsWhat(indexHeap(i)) = 2 * i + 1;

                    swap = valuesHeap(2*i+1);
                    valuesHeap(2*i+1) = valuesHeap(i);
                    valuesHeap(i) = swap;

                    swap = indexHeap(2*i+1);
                    indexHeap(2*i+1) = indexHeap(i);
                    indexHeap(i) = swap;

                    i = 2 * i + 1;
                    continue
                end
            end
            break
        end
    end

    %%


    curIndex = HEAPMinIndex;
    isKnown(curIndex) = true;


    % Mark all neighbours of the newly KNOWN node and mark them as TRIAL
    % (if they are not yet KNOWN).
    newIndices = zeros(numel(HEAPMin), 2*dim);
    for d = 1:dim
        newIndices(:, 2*d-1) = neighborIndices(2+2*(d - 1), curIndex);
        newIndices(:, 2*d) = neighborIndices(1+2*(d - 1), curIndex);
    end
    newIndices(newIndices < 0.5) = [];


    % Update the distance values for all new TRIAL nodes.
    [updatedDistance, updatedOrthogonalS] = updateDistanceWithS(h, 1:dim, ...
        distance, newIndices(:), orthogonalS, neighborIndices);

    old = distance(newIndices(:));
    %distance(newIndices(:)) = updatedDistance;

    for j = 1:numel(newIndices)
        if distance(newIndices(j)) > updatedDistance(j)
            distance(newIndices(j)) = updatedDistance(j);
            orthogonalS(newIndices(j)) = updatedOrthogonalS(j);
        end
    end

    gotUpdated = (abs(updatedDistance - old) > 0);
    changeIndex = newIndices(gotUpdated);
    updateValue = distance(changeIndex);

    %    orthogonalS(changeIndex) = updatedOrthogonalS;

    for j = 1:length(updateValue)

        %%         Update Values in Heap

        NewValue = updateValue(j);
        index = whereIsWhat(changeIndex(j));

        if valuesHeap(index) > uint32(inf) - 1
            kHeap = kHeap + 1;
        end
        valuesHeap(index) = NewValue;

        i = index;

        while abs(i-1) > eps

            temp = uint32(floor(double(i) / 2));
            if (valuesHeap(temp) > valuesHeap(i))

                swap = whereIsWhat(indexHeap(temp));
                whereIsWhat(indexHeap(temp)) = whereIsWhat(indexHeap(i));
                whereIsWhat(indexHeap(i)) = swap;

                swap = valuesHeap(i);
                valuesHeap(i) = valuesHeap(temp);
                valuesHeap(temp) = swap;


                swap = indexHeap(i);
                indexHeap(i) = indexHeap(temp);
                indexHeap(temp) = swap;
                i = temp;
                continue
            end
            break
        end


    end

    if doExit
        break
    end
end

end