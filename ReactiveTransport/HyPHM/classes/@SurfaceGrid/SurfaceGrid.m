%> @brief The class SurfaceGrid represents the triangulation of one-dimensional manifolds @f$\Gamma@f$ on perforated domains.
%>
%> Constructor
%> @code
%>   sg = SurfaceGrid(grid, edgeid)
%> @endcode
%>   where @c grid is an instance of Grid or FoldedGrid and @c edgeid
%>   is the ID (cf. Grid.idE) of some of the edges included in the considered
%>   manifold @f$\Gamma@f$.  The manifold itself is then identified
%>   automatically.

classdef SurfaceGrid

    properties (SetAccess = protected) % only own and derived instances

        %% constants
        %> Underlying (two-dimensional) Grid
        grid
        %> Number of edges of manifold @f$\#E = \#V@f$.
        numE
        %> Number of vertices of manifold @f$\#V = \#E@f$.
        numV

        % Surface edges and surface vertices with original numbers.  Those must
        % subsequently be mapped to 1, 2, 3, ... as being used as indices.
        %> Global edge IDs
        globE
        %> Global vertex IDs
        globV

        %% topology (evaluated when constructor is called)
        %> Vertices to edge @f$[\#V \times \#V]@f$.
        V2E % [#V x #V]
        %> Vertices of edges @f$[\#E \times 2]@f$.
        V0E % [#E x 2 ]

        %% geometry (evaluated when constructor is called)
        %> Coordinates of vertices @f$\vec{x}_V@f$ @f$[\#V \times 2]@f$.
        coordV % [#V x 2]
        %> Length of edges @f$|E|@f$ @f$[\#E \times 1]@f$.
        areaE % [#E x 1 ]    % egde lengths
        %> Barycenters of edges @f$\vec{x}_E^\mathrm{bary}@f$ @f$[\#E \times 2]@f$.
        baryE % [#E x 2 ]    % barycenter of edges
        vecE % [#E x 2 ]    % edge as vector (`tip minus foot')
        %> Unit normals of edges under global orientation @f$\vec{\nu}_E@f$ @f$[\#E \times 2]@f$. The local (triangle-related) normals are @f$\vec{\nu}_{ET}=\sigma_{ET}\vec{\nu}_{E}@f$.
        nuE % [#E x 2 ]    % unit normals on edges
        %> Edge IDs
        idE % [#E x 1]      % edge ids
    end


    methods (Hidden = true) % constructor

        %%% CONSTRUCTOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %> @brief The constructor.
        %> @param varargin @c Grid(grid, edgeid)
        function obj = SurfaceGrid(grid, edgeid)
            assert(isa(grid, 'AbstractGrid'))
            assert(isscalar(edgeid))

            obj.grid = grid;

            % find boundary triangle that contains an edge with the requested ID
            E = find(grid.idE == edgeid);
            E = E(1);
            obj.globE = zeros(0, 1);
            while ~any(obj.globE == E) % as long as E is not contained in the list
                obj.globE(end+1, 1) = E;
                Vright = grid.V0E(E, 2);
                E = find(grid.V0E(:, 1) == Vright); % edges that have Vright as left vertex (also interior edges possible)
                E = E(grid.idE(E) ~= 0); % take the one which is on the boundary (idE ~= 0)
            end
            obj.globV = grid.V0E(obj.globE, 2);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            obj.numV = length(obj.globV);
            obj.numE = length(obj.globE);
            assert(obj.numE == obj.numV);

            % Edges:     |-----1-----|-----2-----|-----3-----|-----4-----|-----5-----|
            % Vertices:  5           1           2           3           4           5
            %%%%%%%%%%%%%%%%%

            %% V2E [#Vx#V] %%
            %%%%%%%%%%%%%%%%%
            % With the numbering above, the matrix V2E [#Vx#V] ("do 2 vertices belong
            % to an edge?") has the following structure:
            %   / 0 2 0 0 1 \
            %   | 2 0 3 0 0 |
            %   | 0 3 0 4 0 |
            %   | 0 0 4 0 5 |
            %   \ 1 0 0 5 0 /
            M = obj.numV;
            obj.V2E = sparse([1, 1:M - 1; M, 2 : M], ...
                [M, 2:M; 1, 1 : M - 1], ...
                [1:M; 1:M]);

            %%%%%%%%%%%%%%%%%

            %% V0E [#Ex#2] %%
            %%%%%%%%%%%%%%%%%
            obj.V0E = [M, 1:M - 1; 1:M]';

            % vertices
            obj.coordV = grid.coordV(obj.globV, :);

            % edge geometry
            obj.baryE = grid.baryE(obj.globE, :);
            obj.vecE = grid.vecE(obj.globE, :);
            obj.idE = grid.idE(obj.globE);
            obj.areaE = grid.areaE(obj.globE);

        end % constructor

    end

end
