%> @file Grid.m The class Grid represents the triangulation of the considered domain @f$\Omega@f$ (or cell domain @f$Y@f$).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @brief The class Grid represents the triangulation of the considered domain @f$\Omega@f$ (or cell domain @f$Y@f$).
%>
%> One possibility to generate an instance of the class Grid is calling
%> @code
%>   g = Grid(coordV, V0T)
%> @endcode
%>   where @c coordV is a list of coordinates of vertices @f$\vec{x}_V@f$
%>   @f$[\#V \times 2]@f$ and @c V0T a list of vertices for each triangle @f$[\#T \times 3]@f$.
%>   The data list describing the mesh may also have the following format:
%>   @f$\mathrm{xyxyxy} = \left[ x_1^i, y_1^i ,x_2^i ,y_2^i ,x_3^i,y_1^i\right]_{i=1,\ldots,\#T}@f$.
%>   Then, simply call
%> @code
%>   g = Grid(xyxyxy)
%> @endcode
%>   A second way is to use the Matlab built-in grid generator @c pdetool and export the
%>   geometric description to the workspace, i.e., @c p, @c e, and @c t. The
%>   constructor can then be called as follows:
%> @code
%>   g = Grid(p, e, t)
%> @endcode
%>   There are several scripts for the generation of triangulations of frequently used domains as rectangles or
%>   perforated squares, which can be found in the directory <tt>HyPHM/gridstuff/generation</tt>.
%>   Some of which are
%> @code
%>   g = domainRectangle(xmin, xmax, ymin, ymax, h)
%>   g = domainCellY(a, b, phi, hmax)
%> @endcode
%>   <b>HyPHM</b> can furthermore directly access <tt>.geo</tt> files of the mesh generator
%>   @ref gmsh or import already generated meshes in the <tt>.mesh</tt> format, e.g.,
%> @code
%>   g = Grid('stdEllipse.geo')  % generates stdEllipse.mesh and extracts data
%> @endcode
%>
%>   Some public methods are  Grid.print (print grid information to
%>   command window),  Grid.ascii (pipe grid information to a file) and
%>   Grid.visualize (plot grid).
%>
%> @image html topo1.jpg
%> @image html topo2.jpg
%>
%> @todo  Deny access for Grid.baryE and Grid.coordV if FoldedGrid (how?).
classdef Grid < AbstractGrid

    methods (Hidden = true) % constructor

        %%% CONSTRUCTOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %> @brief The constructor.
        %> @param varargin @c Grid(coordV, V0T) or @c Grid(initmesh(g)) or @c Grid(p, e, t) or Grid(filename, [scaling])
        function obj = Grid(varargin)
            msg = 'HyPHM: Wrong type of input.  Please see documentation';

            %% GMSH: Grid(initmesh(g)) or Grid(filename, [scaling]) %%%%%%%%%%%%%
            if nargin >= 1 && ischar(varargin{1}) % Grid(filename)
                filename = varargin{1};
                % global mesh size scaling
                if nargin == 2
                    scaling = varargin{2};
                elseif nargin == 1
                    scaling = 1.0;
                else
                    error(msg);
                end
                [~, name, ext] = fileparts(filename);
                assert(strcmp(ext, '.geo') || strcmp(ext, '.mesh'), 'HyPHM: filename has to have extension .geo or .mesh.')

                inputtype = 'gmsh'; % mesh was generated via gmsh
                if strcmp(ext, '.geo') % .geo -> .mesh
                    printline(-1, 'Received .geo file, calling gmsh to convert it to .mesh')
                    geo2mesh(filename, scaling);
                    filename = [name, '.mesh'];
                end

                [coordV, V0T, V0Ebdry, idEbdry] = loadmesh(filename);
                obj.coordV = coordV;
                obj.V0T = V0T;

                %% Grid(xyxyxy) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif nargin == 1 && isnumeric(varargin{1}) % Grid(filename)
                xyxyxy = varargin{1};
                assert(size(xyxyxy, 2) == 6)
                printline(3, 'Initializing Grid.coordV, Grid.V0T (this may take a while)...')
                inputtype = 'manual';
                numT = size(xyxyxy, 1);
                obj.V0T = zeros(numT, 3);
                obj.coordV = zeros(0, 2);
                for kT = 1:numT;
                    for kV = 1:3
                        % compare vertex coordinates
                        marks = obj.coordV(:, 1) == xyxyxy(kT, 2*kV-1) & ... % compare x coordinates
                            obj.coordV(:, 2) == xyxyxy(kT, 2*kV); % compare y coordinates
                        if any(marks) % if we are on a vertex, which is already in our list
                            jV = find(marks);
                            obj.V0T(kT, kV) = jV;
                        else % if we have a new vertex, add the coordinates to our list
                            jV = size(obj.coordV, 1) + 1;
                            obj.coordV(jV, :) = [xyxyxy(kT, 2 * kV - 1), xyxyxy(kT, 2 * kV)];
                            obj.V0T(kT, kV) = jV;
                        end
                    end
                end

                %% Grid(coordV, V0T) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif nargin == 2 % Grid(coordV, V0T)
                assert(isnumeric(varargin{1}), msg)
                assert(isnumeric(varargin{2}), msg)
                inputtype = 'manual';
                obj.coordV = varargin{1};
                obj.V0T = varargin{2};

                %% PDETOOL: Grid(p, e, t) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif nargin == 3 % Grid(initmesh(g)) or Grid(p, e, t)
                for k = 1:3;
                    assert(isnumeric(varargin{k}), msg);
                end
                inputtype = 'pdetool'; % mesh was generated via pdetool from Matlab pde-toolbox
                obj.coordV = varargin{1}';
                obj.V0T = varargin{3}(1:3, :)';
            else
                error(msg)
            end
            %%% now defined: coordV, V0T %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Evaluation of Topology             %
            %   in:  coordV,V0T                 %
            %   out:  sigE0T,V2T,V2E,V0E,T0E,E0T %
            [obj.sigE0T, obj.V2T, obj.V2E, obj.V0E, obj.T0E, ...
                obj.E0T] ...
                = AbstractGrid.evalTopology(obj.coordV, obj.V0T);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Evaluation of local Geometry %
            %   in:  coordV,V0T,E0T,V0E   %
            %   out: coordV0T,baryE0T     %
            [obj.coordV0T, obj.baryE0T] = ...
                AbstractGrid.evalLocGeometry(obj.coordV, obj.V0T, obj.E0T, obj.V0E);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Evaluation of Geometry                  %
            %   in:  coordV,V0T,V0E                   %
            %   out: areaE,baryE,vecE,nuE,areaT,baryT %
            obj = evalGeometry(obj, obj.coordV, obj.V0T, obj.V0E);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Evaluation of Constants   %
            %   in:  V0T,V0E,coordV     %
            %   out: numT,numE,numV     %
            obj.numT = size(obj.V0T, 1);
            obj.numE = size(obj.V0E, 1);
            obj.numV = size(obj.coordV, 1);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Definition of Affine Maps %
            %   in:  coordV0T           %
            %   out: A, b               %

            obj.b = zeros(obj.numT, 2); % [#T x 2]
            for k = 1:2
                obj.b(:, k) = obj.coordV0T(:, 1, k); % coordV0T [#T x 3 x 2]
            end

            obj.A = zeros(obj.numT, 2, 2); % [#T x 2 x 2]
            for k = 1:2
                % first columns of A's
                obj.A(:, k, 1) = obj.coordV0T(:, 2, k) - obj.coordV0T(:, 1, k); % coordV0T [#T x 3 x 2]
                % second columns of A's
                obj.A(:, k, 2) = obj.coordV0T(:, 3, k) - obj.coordV0T(:, 1, k); % coordV0T [#T x 3 x 2]
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%
            % Definition of Edge IDs %
            %   in:  ...             %
            %   out: idE             %
            obj.idE = zeros(obj.numE, 1); % inner edges keep zero
            switch inputtype
                case 'gmsh'
                    for k = 1:size(V0Ebdry, 1)
                        obj.idE(obj.V2E(V0Ebdry(k, 1), V0Ebdry(k, 2))) = idEbdry(k);
                    end
                case 'pdetool'
                    for k = 1:size(varargin{2}, 2)
                        obj.idE(obj.V2E(varargin{2}(1, k), varargin{2}(2, k))) = varargin{2}(5, k);
                    end
                otherwise
                    %printline(-1, 'The edges IDs (Grid.idE) have still to be defined');
            end

        end % constructor

    end


    %%% THE FOLLOWING DOES NOT WORK, SINCE get.baryE %%%
    %%% IS ACCESSED BY CONSTRUCTOR OF PeriodicGrid   %%%
    %   % some get methods
    %   methods
    %     function ret = get.baryE(obj)
    %       assert(~isa(obj, 'PeriodicGrid'), ...
    %         'Using global vertex coordinated is prohibited for periodic Grids.')
    % %       error('foo')
    %       ret = obj.baryE;
    %     end
    %   end

    methods (Hidden = true, Access = private) % called by constructor

        %> Evaluation of geometry of the grid.
        function obj = evalGeometry(obj, coordV, V0T, V0E)
            %        numE0 = size(V0E, 1);
            %        numT0 = size(V0T, 1);
            %       areaE0 = zeros(numE0, 1);
            %       baryE0 = zeros(numE0, 2);
            %       vecE0  = zeros(numE0, 2);
            %       nuE0   = zeros(numE0, 2);
            %       for iEdge = 1 : numE0
            %         vecE0(iEdge, :)  = (coordV(V0E(iEdge, 2), :) - coordV(V0E(iEdge, 1), :))';
            %         areaE0(iEdge)    = norm(vecE0(iEdge, :));
            %         nuE0(iEdge, :)   = vecE0(iEdge, :) * [0,-1; 1,0] / areaE0(iEdge);
            %         baryE0(iEdge, :) = 0.5 * (coordV(V0E(iEdge, 1), :) + coordV(V0E(iEdge,2), :));
            %       end
            vecE0 = (coordV(V0E(:, 2), :) - coordV(V0E(:, 1), :));
            areaE0 = sqrt(vecE0(:, 1).^2+vecE0(:, 2).^2);
            nuE0 = ([vecE0(:, 2), -vecE0(:, 1)] ./ areaE0);
            baryE0 = 0.5 * (coordV(V0E(:, 1), :) + coordV(V0E(:, 2), :));

            % norm(vecE02-vecE0)
            % norm(areaE02-areaE0)
            % norm( nuE0-nuE02)
            % norm(baryE0-baryE02)

            %       areaT0 = zeros(numT0, 1);
            %       baryT0 = zeros(numT0, 2);
            %       for iTriang = 1 : numT0
            %         areaT0(iTriang)    = det([1,1,1; coordV(V0T(iTriang, :), :)'])/2;
            %         baryT0(iTriang, :) = sum(coordV(V0T(iTriang, :), :))/3;
            %       end

            diff1 = coordV(V0T(:, 1), :) - coordV(V0T(:, 2), :);
            diff2 = coordV(V0T(:, 3), :) - coordV(V0T(:, 2), :);
            areaT0 = 0.5 * abs(diff1(:, 1).*diff2(:, 2)-diff1(:, 2).*diff2(:, 1));
            baryT0 = (coordV(V0T(:, 1), :) + coordV(V0T(:, 2), :) + coordV(V0T(:, 3), :)) / 3;

            % norm(areaT0-areaT02)
            % norm(baryT02-baryT0)

            obj.areaE = areaE0;
            obj.baryE = baryE0;
            obj.vecE = vecE0;
            obj.nuE = nuE0;
            obj.areaT = areaT0;
            obj.baryT = baryT0;
        end % evalGeometry

    end % methods

end
