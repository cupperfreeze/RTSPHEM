function [gridNew, levelSetNew] = localMeshRefinementUnFolded(helpGridHyPHM, ...
    levelSet)
% analogue to 'localMeshRefinement' for unfolded grids
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

addedNodes = unique(addedNodes, 'rows');
addedNodes = addedNodes(1:(end-1),:); %delete dummy
deactNodes = unique(deactNodes, 'rows');
deactNodes = deactNodes(2:end); %delete dummy


numaddedNodes = size(addedNodes, 1);
Nodes = [helpGridHyPHM.coordV; addedNodes(:, :)];
DT = delaunayTriangulation(Nodes); %perform Delaunay triangluation of new set of points

GridHyPHM = Grid(DT.Points,DT.ConnectivityList);



nodes = helpGridHyPHM.numV;
levelSet(deactNodes) = eps; %manually deactivated nodes
levelSetNew = [levelSet(:); eps * ones(numaddedNodes, 1)];
gridNew = GridHyPHM;
end