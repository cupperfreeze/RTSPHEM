% Perform local mesh refinement along fluid-solid interface
function [gridNew, levelSetNew] = localMeshRefinement(helpGridHyPHM, ...
    levelSet, varargin)

addedNodes = ones(helpGridHyPHM.numE, 2) * 100; %add nodes if splitting cells is stable
deactNodes = zeros(helpGridHyPHM.numE, 1); %deactivate node if not

indizes = find(logical(levelSet(helpGridHyPHM.V0E(:, 1)) .* levelSet(helpGridHyPHM.V0E(:, 2)) < 100 * eps));
for k = indizes'; %1:helpGridHyPHM.numE
    [deactNodes(k), addedNodes(k, :)] = addpointTripel(helpGridHyPHM.V0E(k, :), ...
        helpGridHyPHM.coordV(helpGridHyPHM.V0E(k, :), :), levelSet(helpGridHyPHM.V0E(k, :)));
end

partitions = round(sqrt(helpGridHyPHM.numV));
addedNodes = unique(addedNodes, 'rows');
deactNodes = unique(deactNodes, 'rows');
deactNodes = deactNodes(2:end); %delete dummy

addedNodes = addedNodes(1:end-1, :);
if any(ismember(varargin, 'bound'))
    sortoutX = logical((addedNodes(:, 1) > -0.5+eps).*(addedNodes(:, 1) < 0.5 - 1 / partitions - eps)); %delete points critical to folding
    sortoutY = logical((addedNodes(:, 2) > -0.5+eps).*(addedNodes(:, 2) < 0.5 - 1 / partitions - eps));
    addedNodes = addedNodes(logical(sortoutX .* sortoutY), :);
end

numaddedNodes = size(addedNodes, 1);
oldNodes = helpGridHyPHM.coordV;

Nodes = [oldNodes; addedNodes(:, :)];
DT = delaunayTriangulation(Nodes);

if any(ismember(varargin, 'Torus'))
    GridHyPHM = FoldedGrid(Grid(DT.Points, DT.ConnectivityList));

    nodes = helpGridHyPHM.numV;
    levelSet(deactNodes) = eps; %manually deactivated nodes
    levelSet = reshape(levelSet, round(sqrt(nodes)), round(sqrt(nodes)));
    levelSet = levelSet(1:round(sqrt(nodes) - 1), 1:round(sqrt(nodes) - 1));
    levelSetNew = [levelSet(:); eps * ones(numaddedNodes, 1)];
    gridNew = GridHyPHM;

elseif any(ismember(varargin, 'Zylinder'))
    GridHyPHM = Grid(DT.Points, DT.ConnectivityList);

    nodes = helpGridHyPHM.numV;
    levelSet(deactNodes) = eps; %manually deactivated nodes
    levelSet = reshape(levelSet, n, m);
    levelSet = levelSet(:, 1:m-1);
    levelSetNew = [levelSet(:); eps * ones(numaddedNodes, 1)];
    gridNew = GridHyPHM;

else
    GridHyPHM = Grid(DT.Points, DT.ConnectivityList);
    levelSet(deactNodes) = eps; %manually deactivated nodes
    levelSetNew = [levelSet(:); eps * ones(numaddedNodes, 1)];
    gridNew = GridHyPHM;

end
end