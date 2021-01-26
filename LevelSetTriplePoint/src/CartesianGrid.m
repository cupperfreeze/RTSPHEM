classdef CartesianGrid < handle
    % Class for cartesian grids on rectangular domains.

    properties (GetAccess = public, SetAccess = immutable)
        dimension; % double scalar
        nodes; % double scalar
        nodesPerDimension; % double array (1 x dim)
        lowerBounds; % double array (1 x dim)
        upperBounds; % double array (1 x dim)
        stepSize; % double array (1 x dim)
    end

    properties (GetAccess = public, SetAccess = public)
        coordinates; % double array (nodes x dimension)

        sub2ind; % double array (size == nodesPerDimension)
        % Physical indices of nodes. This array stores the
        % indices of nodes as seen in the geometry. Even on a
        % folded grid, all nodes are distinct.
        % Example use: coordinates of nodes.

        ind2sub; % double array (nodes x dimension)
        % Saves all subscripts for linear indices.

        triangles; % double array (nodes x dimension+1)
        % Stores indices corresponding to a triangulation of the
        % grid in simplices.

        edges;

        triangleEdges;

        edgeOrientation;

        pull; % short access to neighbors of order 1 and 2
    end


    properties (GetAccess = public, SetAccess = protected)
        isBoundaryEdge;
    end

    methods

        function this = CartesianGrid(dim, bounds, numParts)
            assert(numel(dim) < 2 && numel(dim) > 0);
            assert(abs(dim - round(dim)) < eps);
            assert(dim > 0 && dim < 3, ['Dimension not supported. ', ...
                'Supported dimensions: 1, 2.']);

            this.dimension = dim;

            numParts = numParts(:)';
            this.nodes = prod(numParts+1);

            this.nodesPerDimension = numParts + 1;
            if (this.dimension == uint8(1))
                this.nodesPerDimension = [this.nodesPerDimension, 1];
            end

            this.lowerBounds = bounds(1:2:end);
            this.lowerBounds = this.lowerBounds(:)';

            this.upperBounds = bounds(2:2:end);
            this.upperBounds = this.upperBounds(:)';

            this.stepSize = (this.upperBounds - this.lowerBounds) ...
                ./ numParts;

            this.sub2ind = (1:this.nodes)';
            this.sub2ind = this.reshape(this.sub2ind);

            this.assembleCoordinates();
            this.assembleSubscripts();
            this.createTriangulation();

            %             this.isFoldedGrid = false;
            %             this.index = this.sub2ind;
            %             if ( toFold )
            %                 this.foldGrid();
            %             end

            %speedup in dim=2
            if ~isa(this, 'FoldedCartesianGrid') & abs(dim-2) < eps

                local = cell(2, 2, 2); %direction f/b, depth, dim
                for i = 1:2
                    local{1, 1, i} = this.getForward(1:this.nodes, i);
                    local{2, 1, i} = this.getBackward(1:this.nodes, i);
                    local{1, 2, i} = this.getRelative(1:this.nodes, 2*(1:this.dimension == i));
                    local{2, 2, i} = this.getRelative(1:this.nodes, -2*(1:this.dimension == i));
                end
                this.pull = local;

            end
        end

        % Shifts the nodes with indices in INDEX by the vector SHIFT. If a
        % node with index I has subscripts [X(1), X(2), ..., X(D)], the
        % returned index J denotes the node with subscripts
        % [X(1)+SHIFT(1), X(2)+SHIFT(2), ..., X(D)+SHIFT(D)].
        function retIndex = getRelative(this, index, shift)
            assert(abs(numel(shift) - this.dimension) < eps);

            subsc = this.ind2sub(index, :);
            subsc = subsc + shift(:)';
            hasOverflown = false(numel(index), this.dimension);
            hasUnderflown = false(numel(index), this.dimension);

            for d = 1:this.dimension
                hasOverflown(:, d) = (subsc(:, d) > this.nodesPerDimension(d));
                hasUnderflown(:, d) = (subsc(:, d) < 1);
            end

            retIndex = NaN(numel(index), 1);
            % If this grid is folded, all elements of this array must be TRUE.
            isInterior = all(~hasOverflown & ~hasUnderflown, 2);
            subscCell = mat2cell(subsc(isInterior, :), ...
                sum(isInterior), ones(1, this.dimension));
            if this.dimension ~= 2
                retIndex(isInterior) = diag(this.sub2ind(subscCell{:}));
            else
                retIndex(isInterior) = this.sub2ind(sub2ind(size(this.sub2ind), subscCell{1}, subscCell{2}));
            end
        end

        % Returns the indices of nodes that lie backwards from nodes with
        % index INDEX in the dimension specified by DIM.
        %
        % This result is identical to the function call
        % this.getRelative( INDEX, -1*(1:this.dimension == DIM) ).
        function retIndex = getBackward(this, index, dim)
            assert(abs(numel(dim) - 1) < eps);

            retIndex = this.getRelative(index, -(1:this.dimension == dim));
        end

        % Returns the indices of nodes that lie forwards from nodes with
        % index INDEX in the dimension specified by DIM.
        %
        % This result is identical to the function call
        % this.getRelative( INDEX, 1*(1:this.dimension == DIM) ).
        function retIndex = getForward(this, index, dim)
            assert(abs(numel(dim) - 1) < eps);

            retIndex = this.getRelative(index, (1:this.dimension == dim));
        end

        % Reshapes an array VAL to obtain a structure similar to the
        % geometrical setting of the grid.
        function val = reshape(this, val)
            val = reshape(val, this.nodesPerDimension);
        end

        % Folds the grid so that the grid describes a periodic domain.
        function foldedGrid = foldGrid(this)
            if (abs(this.dimension - 1) < 0.5)
                numParts = this.nodesPerDimension(1) - 1;
            else
                numParts = this.nodesPerDimension - 1;
            end
            bounds = zeros(1, 2*this.dimension);
            bounds(1:2:end) = this.lowerBounds;
            bounds(2:2:end) = this.upperBounds;
            foldedGrid = FoldedCartesianGrid(this.dimension, bounds, numParts);
        end

    end

    methods (Access = private)

        % Creates the matrix of coordinates.
        function assembleCoordinates(this)
            if (abs(this.dimension - 1) < eps)
                this.coordinates = (this.lowerBounds:this.stepSize:this.upperBounds)';
                return;
            end
            intervalCell = cell(1, this.dimension);
            for d = 1:this.dimension
                intervalCell{d} = ...
                    (this.lowerBounds(d):this.stepSize(d):this.upperBounds(d))';
            end
            coordCell = cell(1, this.dimension);
            [coordCell{:}] = ndgrid(intervalCell{:});
            for d = 1:this.dimension
                coordCell{d} = coordCell{d}(:);
            end
            coord = cell2mat(coordCell);
            this.coordinates = coord;
        end

        % Computes the mapping from linear indices to subscripts.
        function assembleSubscripts(this)
            subsc = cell(1, this.dimension);
            % This calls the MATLAB routine ind2sub! (Not the property of this
            % class!)
            [subsc{:}] = ind2sub(this.nodesPerDimension, (1:this.nodes)'); %#ok<CPROP>
            this.ind2sub = cell2mat(subsc);
        end

        % Creates an underlying triangulation of the grid.
        function createTriangulation(this)
            switch (this.dimension)
                case 1
                    this.triangulation1D();
                case 2
                    this.triangulation2D();
                otherwise
                    error('Dimension inconsistency detected.');
            end
        end

        function triangulation1D(this)
            this.triangles = [(1:this.nodes - 1)', (2:this.nodes)'];
            this.edges = (1:this.nodes)';
        end

        % Implementation of a Friedrichs-Keller triangulation. All indices for
        % triangles are stored in counter-clockwise orientation. The first index
        % corresponds to the node on the right angle.
        function triangulation2D(this)
            nX = this.nodesPerDimension(1);
            nY = this.nodesPerDimension(2);
            numSquares = (nX - 1) * (nY - 1);
            sumNumParts = nX + nY - 2;

            triangleIndices = zeros(2*numSquares, 3);
            triIndex = (1:this.nodes)';

            triIndexRed = ...
                triIndex(abs(mod(triIndex, nX)) > 0.5);
            triIndexRed(triIndexRed > nX*(nY - 1)) = [];
            triangleIndices(1:2:end, :) = [triIndexRed, triIndexRed + 1, ...
                triIndexRed + nX];

            triIndexRed = triIndexRed + nX + 1;
            triangleIndices(2:2:end, :) = [triIndexRed, triIndexRed - 1, ...
                triIndexRed - nX];

            this.triangles = triangleIndices;

            numEdges = 3 * numSquares + sumNumParts;
            edgeIndices = zeros(numEdges, 2);
            % The first edges are the edges corresponding to the triangles whose
            % right angle is on the lower left.
            edgeIndices(1:3*numSquares, :) = ...
                [triangleIndices(1:2:end, 1:2); ...
                triangleIndices(1:2:end, 2:3); ...
                triangleIndices(1:2:end, [3, 1])];
            % The next edges are the edges on the right boundary.
            edgeIndices(3*numSquares+1:3*numSquares+nY-1, :) = ...
                nX * [(2:nY)', (1:nY - 1)'];
            % The last edges are the edges on the upper boundary.
            edgeIndices(3*numSquares+nY:end, :) = ...
                nX * (nY - 1) + [(1:nX - 1)', (2:nX)'];

            this.edges = edgeIndices;

            triangleEdgeIndices = zeros(2*numSquares, 3);
            triangleEdgeIndices(1:2:end, :) = ...
                reshape(1:3*numSquares, numSquares, 3);
            triangleEdgeIndices(2:2:end, 1) = ...
                [(nX:numSquares)'; ...
                (3 * numSquares + nY:3 * numSquares + sumNumParts)'];
            triangleEdgeIndices(2:2:end, 2) = ...
                (numSquares + 1:2 * numSquares)';
            vertEdges = zeros(nX, nY-1);
            vertEdges(1:end-1, :) = ...
                reshape(2*numSquares+1:3*numSquares, nX-1, nY-1);
            vertEdges(end, :) = (3 * numSquares + 1:3 * numSquares + nY - 1)';
            vertEdges(1, :) = [];
            triangleEdgeIndices(2:2:end, 3) = vertEdges(:);

            this.triangleEdges = triangleEdgeIndices;

            this.edgeOrientation = kron(ones(numSquares, 3), [1; -1]);

            isBndEdge = false(numEdges, 1);

            edgeCoordX = reshape(this.coordinates(edgeIndices, 1), ...
                size(edgeIndices));
            edgeCoordY = reshape(this.coordinates(edgeIndices, 2), ...
                size(edgeIndices));

            isBndEdge = isBndEdge ...
                | all(abs(edgeCoordX - this.lowerBounds(1)) < eps, 2);
            isBndEdge = isBndEdge ...
                | all(abs(edgeCoordX - this.upperBounds(1)) < eps, 2);
            isBndEdge = isBndEdge ...
                | all(abs(edgeCoordY - this.lowerBounds(2)) < eps, 2);
            isBndEdge = isBndEdge ...
                | all(abs(edgeCoordY - this.upperBounds(2)) < eps, 2);

            this.isBoundaryEdge = isBndEdge;

        end

    end

end
