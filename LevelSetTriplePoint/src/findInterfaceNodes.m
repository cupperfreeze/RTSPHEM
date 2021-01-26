function [isInitial] = findInterfaceNodes(grid, lsf)
% Finds the nodes of a cartesian grid that are adjacent to the interface LSF=0.

assert(isa(grid, 'CartesianGrid'), ...
    'Argument not of class ''CartesianGrid''.');

numNodes = grid.nodes;
isInitial = false(numNodes, 1);
val = NaN(3, numNodes);
val(1, :) = lsf(:)';

for d = 1:grid.dimension

    forwardIndex = grid.pull{1, 1, d};
    val(2, ~isnan(forwardIndex)) = ...
        lsf(forwardIndex(~isnan(forwardIndex)));

    backwardIndex = grid.pull{2, 1, d};
    val(3, ~isnan(backwardIndex)) = ...
        lsf(backwardIndex(~isnan(backwardIndex)));

    isInitial(min(val).*max(val) <= 0) = true;
end


end
