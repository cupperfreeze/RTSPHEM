function [globalStiffness, rightHandSide, isDegreeOfFreedom, ...
    triangleVolumes, triangleSurfaces] = ...
    assembleCellProblem(grid, levelSet)
% Assemble FE system for 2D cell problems and return area of each triangle in fluid (triangleVolumes)
% and length of interface contained in each triangle(triangleSurfaces).

EPS = eps;

assert(isa(grid, 'CartesianGrid'));
isFoldedGrid = isa(grid, 'FoldedCartesianGrid');

numTriangles = size(grid.triangles, 1);
numEdges = size(grid.edges, 1);
numNodes = grid.nodes;
%     maxNumGlobalStiffness = numNodes^2;
maxNumGlobalStiffness = 10 * numNodes;
isIrrelevantNode = true(numNodes, 1);

% TODO: Koeffizient anpassen
rowColIndex = NaN(2*maxNumGlobalStiffness, 2);

matrixEntries = zeros(2*maxNumGlobalStiffness, 1);
rightHandSide = zeros(numNodes, 2);

% Local stiffness matrix of linear FE on triangles.
localStiffness = [2, -1, -1; -1, 1, 0; -1, 0, 1];
referenceVolume = 0.5;

checkedEdge = false(numEdges, 1);
interfaceEdgeCut = NaN(numEdges, 1);
edgeCutPoints = NaN(numEdges, 2);


% In any cut triangle, there must be an edge E with interfaceEdgeCut(E) =
% NaN. This array maps this edge to the node of the triangle, which also is
% a vertex of the triangle that results from the cut with the interface.
% This also doubles as the local index of the edge to the left of a node.
leftEdgeOfVertex = [3; 1; 2];
rightEdgeOfVertex = [2; 3; 1];

interfaceLengths = zeros(numTriangles, 1);
interfaceNormals = NaN(numTriangles, 2);
triangleVolumes = zeros(numTriangles, 1);
triangleSurfaces = zeros(numTriangles, 1);

%%                                                                                  % Handle triangles with trivial interface intersection vectorized
triangleLocalToGlobalIndex = grid.triangles(:, :);
if (isFoldedGrid)
    triangleLocalToGlobalIndex = grid.foldedIndex( ...
        triangleLocalToGlobalIndex);
end
localLevelSet = levelSet(triangleLocalToGlobalIndex);
CompleteSolid = all(localLevelSet > EPS, 2);
CompleteFluid = all(localLevelSet < -EPS, 2);
IndexCompleteFluid = find(CompleteFluid);
CloseToInterface = ~CompleteSolid & ~CompleteFluid;

% Trivial local assembly if triangle completely in solid domain
triangleVolumes(CompleteFluid) = referenceVolume;
isIrrelevantNode(triangleLocalToGlobalIndex(CompleteFluid, :)) = false;

% Fully vectorized local assembly if triangle completely in fluid domain
rowNodes = repmat(triangleLocalToGlobalIndex(IndexCompleteFluid, :), [1, 3]);
colNodes = kron(triangleLocalToGlobalIndex(IndexCompleteFluid, :), [1, 1, 1]);

if (isFoldedGrid)
    xIsFolded = grid.isFoldedNode(grid.triangles(rowNodes));
    rowNodes(xIsFolded) = grid.foldedIndex(rowNodes(xIsFolded));
    yIsFolded = grid.isFoldedNode(grid.triangles(colNodes));
    colNodes(yIsFolded) = grid.foldedIndex(colNodes(yIsFolded));
end

rowColIndex(1:numel(colNodes), 1) = rowNodes(:);
rowColIndex(1:numel(colNodes), 2) = colNodes(:);
matrixEntries(1:numel(colNodes)) = referenceVolume * kron(localStiffness(:), ones(size(colNodes, 1), 1));

indexCounter = numel(colNodes) + 1;


% Treatment of interface triangles
for tri = find(CloseToInterface)'

    triVolume = 0;
    triangleLocalToGlobalIndex = grid.triangles(tri, :);
    if (isFoldedGrid)
        triangleLocalToGlobalIndex = grid.foldedIndex( ...
            triangleLocalToGlobalIndex);
    end
    localLevelSet = levelSet(triangleLocalToGlobalIndex);

    numPositiveEdgeValues = sum(localLevelSet > EPS);
    numNegativeEdgeValues = sum(localLevelSet < -EPS);
    numZeroEdgeValues = 3 - numPositiveEdgeValues - numNegativeEdgeValues;
    %         assert( numZeroEdgeValues < 3 - EPS , ...
    %             'Zero contour has nonzero Lebesgue measure.' );
    if (numZeroEdgeValues > 3 - EPS)
        continue;
    end

    if (numZeroEdgeValues > 2 - EPS)
        % One side of triangle is on the interface.
        % Assemble Neumann boundary condition for this edge.

        nonzeroIndex = find(abs(localLevelSet) > EPS);
        if (sign(localLevelSet(nonzeroIndex)) > EPS)
            % The interior of the triangle lies completely in the positive
            % domain. As such, this triangle will be disregarded (this also
            % assures that the interface boundary condition is evaluated
            % only once.
            continue;
        end
        isIrrelevantNode(triangleLocalToGlobalIndex) = false;
        if (abs(nonzeroIndex - 1) < EPS)
            interfaceLengths(tri) = norm(grid.stepSize, 2);
            normalSign = grid.edgeOrientation(tri, 2);
            interfaceNormals(tri, :) = normalSign ...
                * [grid.stepSize(2), grid.stepSize(1)] ...
                ./ interfaceLengths(tri);
        elseif (abs(nonzeroIndex - 2) < EPS)
            interfaceLengths(tri) = grid.stepSize(2);
            normalSign = -grid.edgeOrientation(tri, 3);
            interfaceNormals(tri, :) = [normalSign, 0];
        else
            interfaceLengths(tri) = grid.stepSize(1);
            normalSign = -grid.edgeOrientation(tri, 1);
            interfaceNormals(tri, :) = [0, normalSign];
        end

        % Midpoint rule for the linear basis functions on the edges lead to
        % coefficients 0.5 (edge adjacent to the corresponding node) or 0
        % (edge opposite to the node).
        localRhs = 0.5 * interfaceLengths(tri) ...
            * repmat(-interfaceNormals(tri, :), 3, 1);
        % The triangle edges are enumerated such that the first vertex of
        % edge I is node I. Thus the opposite edge is edge I+1 (with
        % circular shift).
        %             rightIndex = rightEdgeOfVertex( nonzeroIndex );
        %             rightIndex
        %             localRhs( rightIndex, : ) = 0;
        localRhs(nonzeroIndex, :) = 0;
        rightHandSide(triangleLocalToGlobalIndex, :) = ...
            rightHandSide(triangleLocalToGlobalIndex, :) + localRhs;

        triVolume = referenceVolume;

    elseif (max(localLevelSet) * min(localLevelSet) < -EPS)
        %             tri
        %             max( localLevelSet ) * min( localLevelSet )
        % If at most one vertex is on the zero contour, then the interface
        % doesn't cut the triangle at all or partitions it nontrivially.
        % The latter only happens if there are vertices on both subdomains.
        % If the former is true, then either the standard discretization is
        % used (if negative) or the triangle is disregarded (if positive).
        isIrrelevantNode(triangleLocalToGlobalIndex) = false;

        edges = grid.triangleEdges(tri, :);

        % Calculate cuts of grid edges and the interface.
        for le = 1:3

            e = edges(le);

            if (checkedEdge(e))
                continue;
            end

            checkedEdge(e) = true;

            %phiLeft = levelSet( grid.edges( e, 1 ) );
            %phiRight = levelSet( grid.edges( e, 2 ) );
            phiLeft = levelSet(grid.foldedIndex(grid.edges(e, 1)));
            phiRight = levelSet(grid.foldedIndex(grid.edges(e, 2)));

            % Both nodes lie on the same side of the interface.
            if (phiLeft * phiRight > EPS^2)
                continue;
            end

            % Both nodes lie on the interface.
            if (abs(phiLeft) < EPS && abs(phiRight) < EPS)
                interfaceEdgeCut(e) = Inf;
                continue;
            end

            if (abs(phiLeft - phiRight) < EPS)
                interfaceEdgeCut(e) = Inf;
                continue;
            end

            interfaceEdgeCut(e) = phiLeft / (phiLeft - phiRight);

            assert(interfaceEdgeCut(e) > -EPS & interfaceEdgeCut(e) < 1+EPS);

            edgeCutPoints(e, :) = (1 - interfaceEdgeCut(e)) ...
                * grid.coordinates(grid.edges(e, 1), :) + interfaceEdgeCut(e) ...
                * grid.coordinates(grid.edges(e, 2), :);
        end

        localEdgeCuts = interfaceEdgeCut(edges);
        localCutPoints = edgeCutPoints(edges, :);

        % If the triangle is cut by the interface, there must be at least
        % two edges that are cut. Thus if two or more edges have NaN in
        % their entry of interfaceEdgeCut, there must be an error.
        assert(sum(isnan(localEdgeCuts)) < 1.1);
        if (sum(isnan(localEdgeCuts)) < 0.1)

            % If all edges are cut by the interface, this means that it must
            % cut through a vertex. As such in the current implementation
            % there must be edges with interfaceEdgeCut = 0 or = 1.
            assert(any(abs(localEdgeCuts) < EPS) ...
                && any(abs(localEdgeCuts - 1) < EPS));

            triIndex = find(abs(localLevelSet - min(localLevelSet)) < eps);
        else
            triIndex = leftEdgeOfVertex(isnan(localEdgeCuts));
        end

        leftIndex = leftEdgeOfVertex(triIndex);
        rightIndex = rightEdgeOfVertex(triIndex);

        volumeTrianglePart = 0.5 * ...
            (0.5 - grid.edgeOrientation(tri, leftIndex) * ...
            (localEdgeCuts(leftIndex) - 0.5)) * ...
            (0.5 + grid.edgeOrientation(tri, triIndex) * ...
            (localEdgeCuts(triIndex) - 0.5));

        if (localLevelSet(triIndex) < EPS)
            triVolume = volumeTrianglePart;
        else
            triVolume = 0.5 - volumeTrianglePart;
        end

        interfaceLengths(tri) = norm(localCutPoints(triIndex, :) ...
            -localCutPoints(leftIndex, :));
        if (abs(interfaceLengths(tri)) < EPS)
            break;
        end
        interfaceNormals(tri, :) = (localCutPoints(leftIndex, :) ...
            -localCutPoints(triIndex, :)) * [0, -1; 1, 0];
        interfaceNormals(tri, :) = interfaceNormals(tri, :) ...
            ./ interfaceLengths(tri);
        if (localLevelSet(triIndex) > EPS)
            interfaceNormals(tri, :) = -interfaceNormals(tri, :);
        end

        localRhs = repmat(-interfaceNormals(tri, :), 3, 1);
        localRhs(triIndex, :) = localRhs(triIndex, :) ...
            * 0.5 * interfaceLengths(tri) ...
            * ((0.5 - grid.edgeOrientation(tri, leftIndex) * ...
            (localEdgeCuts(triIndex) - 0.5)) ...
            +(0.5 + grid.edgeOrientation(tri, triIndex) * ...
            (localEdgeCuts(leftIndex) - 0.5)));
        localRhs(leftIndex, :) = localRhs(leftIndex, :) ...
            * 0.5 * interfaceLengths(tri) ...
            * (0.5 - grid.edgeOrientation(tri, leftIndex) * ...
            (localEdgeCuts(leftIndex) - 0.5));
        localRhs(rightIndex, :) = localRhs(rightIndex, :) ...
            * 0.5 * interfaceLengths(tri) ...
            * (0.5 + grid.edgeOrientation(tri, triIndex) * ...
            (localEdgeCuts(triIndex) - 0.5));
        rightHandSide(triangleLocalToGlobalIndex, :) = ...
            rightHandSide(triangleLocalToGlobalIndex, :) + localRhs;

    elseif (max(localLevelSet) < EPS)
        isIrrelevantNode(triangleLocalToGlobalIndex) = false;
        triVolume = referenceVolume;
    end

    % [x,y] = meshgrid( grid.triangles( tri, : ) );
    x = repmat(grid.triangles(tri, :), [3, 1]);
    y = x';

    if (isFoldedGrid)
        %  [ xIsFolded, yIsFolded ] = ...
        %      meshgrid( grid.isFoldedNode( grid.triangles(tri,:) ) );
        xIsFolded = repmat(grid.isFoldedNode(grid.triangles(tri, :))', [3, 1]);
        yIsFolded = xIsFolded';
        x(xIsFolded) = grid.foldedIndex(x(xIsFolded));
        y(yIsFolded) = grid.foldedIndex(y(yIsFolded));
    end

    for i = 1:numel(x)
        % index = find( all( bsxfun( @eq, rowColIndex, [ x(i) y(i) ] ), 2 ) );

        %             if ( numel( index ) < 1 - EPS )
        %                 rowColIndex( indexCounter, : ) = [ x(i) y(i) ];
        %                 index = indexCounter;
        %                 indexCounter = indexCounter + 1;
        %             else
        %                 assert( numel( index ) < 2 - EPS );
        %             end

        % rowColIndex( indexCounter, : ) = [ x(i) y(i) ];
        rowColIndex(indexCounter, 1) = x(i);
        rowColIndex(indexCounter, 2) = y(i);
        matrixEntries(indexCounter) = triVolume * localStiffness(i);
        %             index = indexCounter;
        indexCounter = indexCounter + 1;

        %             matrixEntries( index ) = matrixEntries( index ) ...
        %                 + triVolume * localStiffness(i);
    end


    triangleVolumes(tri) = triVolume;

end

triangleVolumes = triangleVolumes .* prod(grid.stepSize);


triangleSurfaces = interfaceLengths;

%     rowColIndex( indexCounter:end, : ) = [];
%     matrixEntries( indexCounter:end ) = [];


toDelete = isIrrelevantNode;
if (isFoldedGrid)
    toDelete = toDelete | grid.isFoldedNode;
end

% global assembly
globalStiffness = sparse(rowColIndex(1:indexCounter - 1, 1), rowColIndex(1:indexCounter - 1, 2), ...
    matrixEntries(1:indexCounter - 1), numNodes, numNodes);

globalStiffness(toDelete, :) = [];
globalStiffness(:, toDelete) = [];
isDegreeOfFreedom = ~toDelete;

if (isFoldedGrid)
    globalStiffness(end+1, :) = 1;
end

end
