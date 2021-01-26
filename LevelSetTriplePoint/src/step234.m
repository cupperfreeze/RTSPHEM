function [dV, S] = step234(distFunctions, grid, Xi, interfaceVelocities, restrictDist)

numNodes = grid.nodes;

isInitial = false(numNodes, 1);
val = NaN(3, numNodes);
val(1, :) = Xi';
initialS = inf(size(Xi, 1), grid.dimension);

%Label all Interface points by checking wheater forward/backward steps in
%any dimension change Xi
for d = 1:grid.dimension

    forwardIndex = grid.pull{1, 1, d};
    val(2, ~isnan(forwardIndex)) = ...
        Xi(forwardIndex(~isnan(forwardIndex)));

    backwardIndex = grid.pull{2, 1, d};
    val(3, ~isnan(backwardIndex)) = ...
        Xi(backwardIndex(~isnan(backwardIndex)));

    %when Xi changes in forward step
    changeIndex = val(1, :) ~= val(2, :) & ~isnan(val(1, :)) & ~isnan(val(2, :));
    isInitial(changeIndex) = true;
    initialS(changeIndex, d) = interfaceVelocities(sub2ind(size(interfaceVelocities), val(1, changeIndex), val(2, changeIndex)));

    %when Xi changes in backward step
    changeIndex = val(1, :) ~= val(3, :) & ~isnan(val(1, :)) & ~isnan(val(3, :));
    isInitial(changeIndex) = true;
    initialS(changeIndex, d) = interfaceVelocities(sub2ind(size(interfaceVelocities), val(1, changeIndex), val(3, changeIndex)));
end


%Assings lowest distances at interface points
[dMinValue, dMinIndex] = min(distFunctions, [], 2);
dVInitial = inf(size(isInitial));
Initials = find(isInitial);
for i = 1:length(Initials)
    InitialIndex = Initials(i);
    dVInitial(InitialIndex) = min(abs(dMinValue(InitialIndex) - distFunctions(InitialIndex, mod((dMinIndex(InitialIndex)+0), 3) + 1)), ...
        abs(dMinValue(InitialIndex) - distFunctions(InitialIndex, mod((dMinIndex(InitialIndex)+1), 3) + 1))) / 2;

    if abs(abs(dMinValue(InitialIndex) - distFunctions(InitialIndex, mod((dMinIndex(InitialIndex)+0), 3) + 1))- ...
            abs(dMinValue(InitialIndex) - distFunctions(InitialIndex, mod((dMinIndex(InitialIndex)+1), 3) + 1))) < 10 * eps;

        initialS(InitialIndex) = sign(initialS(InitialIndex)) * min(abs(interfaceVelocities(dMinIndex(InitialIndex), mod((dMinIndex(InitialIndex)+0), 3) + 1)), ...
            abs(interfaceVelocities(dMinIndex(InitialIndex), mod((dMinIndex(InitialIndex)+1), 3) + 1)));
    end


end

exceptions = find(~isinf(initialS(:, 1)) & ~isinf(initialS(:, 2)) & ...
    initialS(:, 1) ~= initialS(:, 2));
for i = 1:length(exceptions)
    neighbors = [grid.pull{1, 1, 1}(exceptions(i)), grid.pull{2, 1, 1}(exceptions(i)), ...
        grid.pull{1, 1, 2}(exceptions(i)), grid.pull{2, 1, 2}(exceptions(i))];
    if sum(isnan(neighbors)) == 0
        CrossInterfaceNeighbors = neighbors(Xi(exceptions(i)) ~= Xi(neighbors));
        [~, index] = max(dVInitial(CrossInterfaceNeighbors)-eps*abs(initialS(CrossInterfaceNeighbors))');
        index = CrossInterfaceNeighbors(index);
        initialS(exceptions(i), 1) = interfaceVelocities(Xi(exceptions(i)), Xi(index));
    end
end

initialS(isinf(initialS(:, 1)), 1) = initialS(isinf(initialS(:, 1)), 2);
[dV, S] = reinitializeLevelSetWithS(grid, dVInitial, restrictDist, initialS(:, 1)');
S(isinf(S)) = 0;

end