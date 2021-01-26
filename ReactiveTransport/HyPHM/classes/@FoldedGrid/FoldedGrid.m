%> @file FoldedGrid.m A rectangular grid (possibly with insections) with semi or full periodic boundaries.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%>  A rectangular grid (possibly with insections) with semi or full periodic boundaries.
classdef FoldedGrid < AbstractGrid

    properties (SetAccess = private, Hidden)
        %> Flag. If set, the rectangular bounded AbstractGrid is periodic/folded in vertical direction only.
        isCylinder

        %> Associated rectangular-bounded Grid.
        rectGrid
        %> Vertex numbers related to north boundary on rectangular grid @f$ [1\times E_\mathrm{north}] @f$ .
        northV
        eastV % these properties are initialized via getOuterBoundaryVerts().
        southV % northV and southV are ordered wrt x1-direction,
        westV % eastV and westV are ordered wrt x2-direction.

        northE
        eastE
        southE
        westE

        %> @f$[\#V^\mathrm{folded}\times 1]@f$ mapping which descripes the mapping
        %> of the vertices on both grids. The kth entry of mapV @f$[\#V^\mathrm{rectangular}\times 1]@f$
        %> contains the number of the new/periodic vertex (wrt the new numbering)
        %> to which the old/rectangular vertex k (wrt the old/rectangular numbering) is mapped to.
        %> This mapping is <i>not</i> an injection.
        mapV

        %> Analogeous to mapV.  mapT is not needed, since it is the identity.
        mapE

        % mapT is not needed, since it is the identity

    end % private properties

    methods (Hidden = true) % constructor

        %> @brief Constructor which requires a yet complete Grid.
        %> @param grid Rectangular bounded instance of Grid.
        %> @param varargin Void for torus, one for (vertically) periodic cylinder.
        function this = FoldedGrid(grid, varargin)
            assert(isa(grid, 'Grid'))
            if isempty(varargin)
                this.isCylinder = false;
                %printline(3, 'Grid is folded to a torus.')
            else
                this.isCylinder = true;
                %printline(3, 'Grid is folded to a cylinder.')
            end


            this.rectGrid = grid; % storing the original rectangular-bounded grid

            this.searchOuterBoundary; % init of northV, eastV, southV, westV

            % evaluation of basic topology/geometry
            switch this.isCylinder
                case false
                    [coordV] = this.fold2torus; % init of coordV, mapV
                case true
                    [coordV] = this.fold2cylinder;
            end

            % fetching geometry data with local indexing from rectGrid
            % no index mapping has to be performed since i->i for Triangles
            this.coordV0T = this.rectGrid.coordV0T;
            this.baryE0T = this.rectGrid.baryE0T;
            this.A = this.rectGrid.A;
            this.b = this.rectGrid.b;

            % building new triangle-vertex-table
            this.V0T = this.mapV(this.rectGrid.V0T);

            % setting known constants
            this.numV = size(this.coordV, 1);
            this.numT = this.rectGrid.numT;

            % evaluation of topology
            [this.sigE0T, this.V2T, this.V2E, this.V0E, this.T0E, this.E0T] ...
                = AbstractGrid.evalTopology(coordV, this.V0T);

            % setting known constants
            this.numE = size(this.V0E, 1);

            %%% BUILDUNG UP THE MAPPING MAP_EDGES
            this.mapE = zeros(this.rectGrid.numE, 1);
            %       for iRectEdge = 1 : this.rectGrid.numE
            %         [flag, idx] = ismember(this.mapV(this.rectGrid.V0E(iRectEdge, 1:2))', this.V0E(:, 1:2), 'rows');
            %         if flag == 0 % if entry is not found, the edge has a different orientation
            %           [flag, idx] = ismember(this.mapV(this.rectGrid.V0E(iRectEdge, 1:2))', this.V0E(:, [2, 1]), 'rows');  %#ok<ASGLU>
            %         end
            %         this.mapE(iRectEdge) = idx;
            %       end

            [~ idx] = ismember(this.mapV(this.rectGrid.V0E(:, 1:2)), this.V0E(:, 1:2), 'rows'); %improved SG
            [~ idx2] = ismember(this.mapV(this.rectGrid.V0E(:, 1:2)), this.V0E(:, [2, 1]), 'rows');
            idx(idx == 0) = idx2(idx == 0);
            this.mapE(:) = idx;

            %%% EVALUATION OF THE GEOMETRY (w/o coordV)
            this = this.evalPeriodicGeometry;

            %%% building edge ID list
            this.idE = zeros(this.numE, 1);
            this.idE(this.mapE) = this.rectGrid.idE;
            switch this.isCylinder
                case false
                    this.idE(this.mapE(this.northE)) = 0; % delete north and south ids
                    this.idE(this.mapE(this.eastE)) = 0; % delete east and west ids
                case true
                    this.idE(this.mapE(this.northE)) = 0; % delete north and south ids
            end
        end % constructor

    end

end
