function [gridNew, levelSetNew] = localMeshRefinement(helpGridHyPHM, ...
    levelSet)

% adapt the Grid such that edges are aligned with the linear approximation
% of 0-level-set. As the new grid may contain a different number of
% vertices an updates level-set function is returned as well.

indizes = find(logical(levelSet(helpGridHyPHM.V0E(:, 1)) .* ... %edges that are possibly cut by interface
    levelSet(helpGridHyPHM.V0E(:, 2)) < 100 * eps));
addedNodes = inf(numel(indizes)+1, 2); %add nodes if splitting cells is stable
deactNodes = zeros(numel(indizes)+1, 1); %deactivate node if not

indizes = find(logical(levelSet(helpGridHyPHM.V0E(:, 1)) .* levelSet(helpGridHyPHM.V0E(:, 2)) < 100 * eps));
for k = 1:numel(indizes); %1:helpGridHyPHM.numE
    locidx = indizes(k);
    [deactNodes(k), addedNodes(k, :)] = addpoint(helpGridHyPHM.V0E(locidx, :), ...
        helpGridHyPHM.coordV(helpGridHyPHM.V0E(locidx, :), :), levelSet(helpGridHyPHM.V0E(locidx, :)));

end

partitions = round(sqrt(helpGridHyPHM.numV));
addedNodes = unique(addedNodes, 'rows');
deactNodes = unique(deactNodes, 'rows');
deactNodes = deactNodes(2:end); %delete dummy

sortoutX = logical((addedNodes(:, 1) > -0.5+eps).*(addedNodes(:, 1) < 0.5 - 1 / partitions - eps)); %delete points critical to folding
sortoutY = logical((addedNodes(:, 2) > -0.5+eps).*(addedNodes(:, 2) < 0.5 - 1 / partitions - eps));
addedNodes = addedNodes(logical(sortoutX .* sortoutY), :);

numaddedNodes = size(addedNodes, 1);
Nodes = [helpGridHyPHM.coordV; addedNodes(:, :)];
DT = delaunayTriangulation(Nodes); %perform Delaunay triangluation of new set of points

GridHyPHM = FoldedGrid(Grid(DT.Points, DT.ConnectivityList));
%GridHyPHM = Grid(DT.Points,DT.ConnectivityList);
%GridHyPHM.visualize('numV');


nodes = helpGridHyPHM.numV;
levelSet(deactNodes) = eps; %manually deactivated nodes
levelSet = reshape(levelSet, round(sqrt(nodes)), round(sqrt(nodes)));
levelSet = levelSet(1:round(sqrt(nodes) - 1), 1:round(sqrt(nodes) - 1));
levelSetNew = [levelSet(:); eps * ones(numaddedNodes, 1)];
gridNew = GridHyPHM;
end