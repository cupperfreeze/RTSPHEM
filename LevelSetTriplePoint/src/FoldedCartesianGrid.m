classdef FoldedCartesianGrid < CartesianGrid
    % Class for cartesian grids on rectangular domains with periodic boundaries.


    properties (GetAccess = public, SetAccess = immutable)
        foldedIndex; % double array (nodes x 1)
        % Logical indices of nodes. If the grid is folded, some
        % nodes are regarded as identical and are thus assigned
        % the same index.
        % Example use: evaluation of values defined on nodes.
        isFoldedNode; % logical array (nodes x 1)
        realNodes;

        %pull;
    end

    methods

        function this = FoldedCartesianGrid(dim, bounds, numParts)
            this@CartesianGrid(dim, bounds, numParts);

            this.isBoundaryEdge = false(size(this.isBoundaryEdge));
            this.foldedIndex = this.sub2ind;
            this.isFoldedNode = false(this.nodes, 1);
            semicolons = num2cell(repmat(':', 1, this.dimension));
            for d = 1:dim
                this.foldedIndex(semicolons{1:d - 1}, this.nodesPerDimension(d), ...
                    semicolons{d + 1:end}) = ...
                    this.foldedIndex(semicolons{1:d - 1}, 1, ...
                    semicolons{d + 1:end});
                this.isFoldedNode(this.sub2ind(semicolons{1:d - 1}, this.nodesPerDimension(d), ...
                    semicolons{d + 1:end})) = true;
            end

            this.realNodes = sum(~this.isFoldedNode);
            %speedup in dim=2
            if abs(dim-2) < eps
                local = cell(2, 2, 2); %direction f/b, depth, dim
                for i = 1:2
                    local{1, 1, i} = this.getForward(1:this.nodes, i);
                    local{2, 1, i} = this.getBackward(1:this.nodes, i);
                    local{1, 2, i} = this.getForward(local{1, 1, i}, i);
                    local{2, 2, i} = this.getBackward(local{2, 1, i}, i);
                end
                this.pull = local;
            end
        end

        % Shifts the nodes with indices in INDEX by the vector SHIFT. If a
        % node with index I has subscripts [X(1), X(2), ..., X(D)], the
        % returned index J denotes the node with subscripts
        % [X(1)+SHIFT(1), X(2)+SHIFT(2), ..., X(D)+SHIFT(D)].
        function retIndex = getRelative(this, index, shift)

            global EPS;
            EPS = eps;
            assert(abs(numel(shift) - this.dimension) < EPS);

            subsc = this.ind2sub(index, :);
            % The maximal value for a shift in dimension D can be
            % this.nodesPerDimension(D)-1. Greater shifts must lie outside of
            % the domain for any node on non-folded grids. On a folded grid,
            % a greater shift circles the dimension more than once.
            if (this.dimension > 1 + EPS)
                shift = mod(shift, this.nodesPerDimension-1);
            else
                shift = mod(shift, this.nodesPerDimension(1)-1);
            end
            subsc = subsc + shift(:)';
            hasOverflown = false(numel(index), this.dimension);
            hasUnderflown = false(numel(index), this.dimension);

            for d = 1:this.dimension
                hasOverflown(:, d) = (subsc(:, d) > this.nodesPerDimension(d));
                hasUnderflown(:, d) = (subsc(:, d) < 1);
                subsc(hasOverflown(:, d), d) = ...
                    subsc(hasOverflown(:, d), d) ...
                    -this.nodesPerDimension(d) + 1;
                subsc(hasUnderflown(:, d), d) = ...
                    subsc(hasUnderflown(:, d), d) ...
                    +this.nodesPerDimension(d) - 1;
                hasOverflown(:, d) = false;
                hasUnderflown(:, d) = false;
            end

            % If this grid is folded, all elements of this array must be TRUE.
            %isInterior = all( ~hasOverflown & ~hasUnderflown, 2 );
            if (this.dimension > 1 + EPS)
                %                 subscCell = mat2cell( subsc, numel(index), ...
                %                     ones( 1, this.dimension ) );
                %                 retIndex = diag( this.sub2ind( subscCell{:} ) );

                retIndex = subsc(:, 1) + this.nodesPerDimension(1) * (subsc(:, 2) - 1);

                %               General case:
                %                 retIndex = subsc( :, this.dimension ) - 1;
                %                 for d = this.dimension - 1 : 1
                %                     if ( d < 1 )
                %                         break;
                %                     end
                %                     retIndex = this.nodesPerDimension( d ) * retIndex + ...
                %                         subsc( :, d );
                %                 end
            else
                retIndex = subsc;
            end
            retIndex = this.foldedIndex(retIndex);
        end

        function val = synchronizeValues(this, val)
            val(this.isFoldedNode) = ...
                val(this.foldedIndex(this.isFoldedNode));
        end

    end

end
