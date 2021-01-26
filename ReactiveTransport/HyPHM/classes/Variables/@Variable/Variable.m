%> @file Variable.m Class for unknowns, coefficients, and parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @brief  Class for unknowns, coefficients, and parameters.
%>   Each problem contains coefficient functions and unknowns living in
%>   different discrete function spaces.  Both have to be defined by instances
%>   of the class Variable.  For example, we want to define a velocity
%>   @f$\vec{v}@f$ and choose @f$\mathbb{P}_2(\mathcal{T})^2@f$ as function space, ie each
%>   component of @f$\vec{v}@f$ is a triangle-wise quadratic polynomial on the
%>   Grid:
%>   @code
%>     v = Variable(g, st, 'velocity', 'P2P2');
%>     p = Variable(g, st, 'pressure', 'P1');
%>     D = Variable(g, st, '',         'P0P0P0P0'); % the name may also be empty
%>   @endcode
%>   where @c g and @c st are instances of Grid and Stepper, respectively.
%>   The third argument contains the name to be used when visualized via
%>   v.visualize.  The fourth arguments are the tags for the spaces
%>   @f$\mathbb{P}_2(\mathcal{T})^2@f$, @f$\mathbb{P}_1(\mathcal{T})@f$ and @f$\mathbb{P}_0(\mathcal{T})^{2,2}@f$, see Variable.
%>
%> @b Remark:
%>
%> - The user should feed instances of Variable <i>only</i> by Variable.setdata.
%> - Variable.checkdata, Variable.setfh2data (opt), Variable.setfh3data (opt), Variable.visualize (opt) have to
%>   be implemented for each Variable.type, ie for each discrete function
%>   space.
%> - Variable.setfh3data is usually faster than Variable.setfh2data.
%>
%>
%> @warning
%> If the data contains @c inf or @c NaN the visualization tool (like
%> Paraview) will likely crash!
%>
%>
%> @todo Tutorial with visualization of a spiral.

classdef Variable < hgsetget & handle

    properties (SetAccess = protected) % only own and derived instances
        %> Underlying Grid.
        grid
        grids %in case grid varies
        %> Underlying Stepper.
        stepper
        %> @brief Name of the discrete function space, see directory <code>classes/+Variables/</code>.
        %> Some are
        %>
        %> <table>
        %> <tr>
        %> <th align="left" width="100"><b>Space</b></th>
        %> <th align="left" width="100"><b>Tag</b></th>
        %> <th align="left" width="150"><b>Dim of Data</b></th>
        %> <th align="left" width="100"><b>Eucl Dim</b></th>
        %> <th align="left" width="400"><b>Description</b></th>
        %> <th align="left" width="200"><b>DOFs</b></th>
        %> </tr>
        %> <tr>
        %>  <td>@f$\mathbb{P}_0(\Omega)@f$</td>
        %>  <td>@c C</td>
        %>  <td>@f$[1]@f$</td>
        %>  <td>scalar</td>
        %>  <td>globally constant</td>
        %>  <td>one global</td>
        %> </tr>
        %> <tr>
        %>  <td>@f$\mathbb{P}_0(\mathcal{T})@f$</td>
        %>  <td><code>P0</code></td>
        %>  <td>@f$[\#T\times 1]@f$</td>
        %>  <td>scalar</td>
        %>  <td>constant on each triangle (i.e. <i>discontinuous</i>)</td>
        %>  <td>triangles</td>
        %> </tr>
        %> <tr>
        %>  <td>@f$\mathbb{P}_0(\mathcal{E})@f$</td>
        %>  <td><code>P0E</code></td>
        %>  <td>@f$[\#E\times 1]@f$</td>
        %>  <td>scalar</td>
        %>  <td>constant on each edge (i.e. <i>discontinuous</i>)</td>
        %>  <td>edges</td>
        %> </tr>
        %> <tr>
        %>  <td>@f$\mathbb{P}_1^\mathrm{c}(\mathcal{T})@f$</td>
        %>  <td>@c P1</td>
        %>  <td>@f$[\#V\times 1]@f$</td>
        %>  <td>scalar</td>
        %>  <td>1st order polynomial on each triangle, globally <i>continuous</i></td>
        %>  <td>vertices (shared)</td>
        %> </tr>
        %> <tr>
        %>  <td>@f$\mathbb{P}_1(\mathcal{T})@f$</td>
        %>  <td>@c P1disc</td>
        %>  <td>@f$[3\#T\times 1]@f$</td>
        %>  <td>scalar</td>
        %>  <td>1st order polynomial on each triangle, globally <i>discontinous</i></td>
        %>  <td>vertices (not shared)</td>
        %> </tr>
        %> <tr>
        %>  <td>@f$\mathbb{CR}(\mathcal{T})@f$</td>
        %>  <td>@c CR</td>
        %>  <td>@f$[\#E\times 1]@f$</td>
        %>  <td>scalar</td>
        %>  <td>Crouzeix&ndash;Raviart space; 1st order polynomial on each triangle, continuous in each edge barycenter but globally <i>discontinuous</i>.
        %>      @f[\mathbb{CR}(\mathcal{T}) := \Big\{v_h:\Omega\rightarrow\mathbb{R};\, \forall T\in\mathcal{T}_h,\, v_h\big|_T\in\mathbb{P}_1(T);\,\forall E\in\mathcal{E}_\Omega,\, \int_E [v_h] = 0\Big\} @f] </td>
        %>  <td>edges (shared)</td>
        %> </tr>
        %> <tr>
        %>  <td>@f$\mathbb{P}_2^\mathrm{c}(\mathcal{T})@f$</td>
        %>  <td>@c P2</td>
        %>  <td>@f$[\#V+\#E\times 1]@f$</td>
        %>  <td>scalar</td>
        %>  <td>2st order polynomial on each triangle, globally <i>continuous</i>.  The DOFs are numbered the following way:
        %>      @f[\Big[\underbrace{\sigma_1,\ldots,\sigma_{\#V}}_{\text{DOFs on vertices}},\underbrace{\sigma_{\#V+1},\ldots,\sigma_{\#V+\#E}}_{\text{DOFs on edges}}\Big]@f]</td>
        %>  <td>vertices and edges (shared)</td>
        %> </tr>
        %> <tr>
        %>  <td>@f$\mathbb{RT}_0(\mathcal{T})@f$</td>
        %>  <td>@c RT0</td>
        %>  <td>@f$[\#E\times 1]@f$</td>
        %>  <td>vector</td>
        %>  <td>global Raviart&ndash;Thomas space; "halve"-linear in each component; continuous normal trace over edges</td>
        %>  <td>edges (shared)</td>
        %> </tr>
        %> <tr>
        %>  <td>@f$\mathbb{P}_0(\mathcal{T})^2@f$</td>
        %>  <td>@c P0P0</td>
        %>  <td>@f$[\#T\times 2]@f$</td>
        %>  <td>vector</td>
        %>  <td>component-wise in @f$\mathbb{P}_0(\mathcal{T})@f$</td>
        %>  <td>triangles</td>
        %> </tr>
        %> <tr>
        %>  <td>@f$\mathbb{P}_2^\mathrm{c}(\mathcal{T})^2@f$</td>
        %>  <td>@c P2P2</td>
        %>  <td>@f$[\#V+\#E\times 2]@f$</td>
        %>  <td>vector</td>
        %>  <td>component-wise in
        %>   @f$\mathbb{P}_2^\mathrm{c}(\mathcal{T})@f$ (see @f$\mathbb{P}_2^\mathrm{c}(\mathcal{T})@f$ for details)</td>
        %>  <td>vertices and edges (shared)</td>
        %> </tr>
        %> <tr>
        %>  <td>@f$\mathbb{P}_1Bubble (\mathcal{T})^2@f$</td>
        %>  <td>@c P1BubbleP1Bubble</td>
        %>  <td>@f$[\#V+\#T\times 2]@f$</td>
        %>  <td>vector</td>
        %>  <td>component-wise in
        %>   @f$\mathbb{P}_1Bubble^\mathrm{c}(\mathcal{T})@f$ </td>
        %>  <td>vertices (shared) and triangles </td>
        %> </tr>
        %> <tr>
        %>  <td>@f$\mathbb{P}(\Omega)^{2,2}@f$</td>
        %>  <td>@c CCCC</td>
        %>  <td>@f$[4]@f$</td>
        %>  <td>tensor</td>
        %>  <td>component-wise globally constant,   one DOF on @f$\Omega@f$ for each component</td>
        %>  <td>one global</td>
        %> </tr>
        %> <tr>
        %>  <td>@f$\mathbb{P}_0(\mathcal{T})^{2,2}@f$</td>
        %>  <td>@c P0P0P0P0</td>
        %>  <td>@f$[\#T\times 4]@f$</td>
        %>  <td>tensor</td>
        %>  <td>component-wise in @f$\mathbb{P}_0(\mathcal{T})@f$</td>
        %>  <td>triangles</td>
        %> </tr>
        %> </table>
        %>
        %> @b Example:
        %> Consider the function @f$\vec{x}\mapsto \cos(2\pi x)\cos(2\pi y)@f$
        %> on @f$[0.25,0.75]^2@f$.  This function, mapped to three discrete spaces, looks like:
        %>
        %> <table>
        %> <tr>
        %>   <td align="center" bgcolor="#FFFFFF"> @image html images/variable-P1.png </td>
        %>   <td align="center" bgcolor="#FFFFFF"> @image html images/variable-CR.png </td>
        %>   <td align="center" bgcolor="#FFFFFF"> @image html images/variable-P0.png </td>
        %> </tr>
        %> <tr>
        %>  <td align="center"> @f$\mathbb{P}_1^\mathrm{c}(\mathcal{T})@f$ </td>
        %>  <td align="center"> @f$\mathbb{CR}(\mathcal{T})@f$ </td>
        %>  <td align="center"> @f$\mathbb{P}_0(\mathcal{T})@f$ </td>
        %> </tr>
        %> </table>
        %>
        %> The underlying code was
        %> @code
        %> g = domainRectangle(.25, .75, .25, .75, .1);
        %> st = Stepper(0:1);
        %> for kS = {'P1', 'CR', 'P0'}
        %>   v = Variable(g, st, ['v-' kS{1}], kS{1});
        %>   v.setdata(@(t, X, Y) cos(2*pi*X).*sin(2*pi*Y));
        %>   v.visualize
        %> end
        %> @endcode
        type
    end

    properties % public access
        %> Description [string], exclusively used for vtk visualization (cf. <code>Variable.visualize</code>).
        name
    end

    properties (Access = private, Hidden)
        %> This is where the internal data is stored, index shift! (1 is initial)
        %> cell array with length equal to numsteps + 1, where one entry is
        content
    end

    %% constructor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        %> The constructor.
        %> @param grid Instance of AbstractGrid
        %> @param stepper Instance of Stepper
        %> @param name Description of the variable (string, optional), if void [] a random name is used.
        %> @param type Discrete function space, see above
        function this = Variable(grid, stepper, name, type)
            msg = ['HyPHM: Variable(Grid, Stepper, name, type).  ', ...
                'See documentation (>>doc Variable) for further information'];
            assert(nargin == 4, msg);
            if ~strcmp(type, 'C') && ~strcmp(type, 'CCCC') % if variable is a constant it is independent from the domain
                this.grid = grid; % verified in set function
            end
            this.grids = cell(stepper.numsteps, 1);
            this.grids{1} = grid;
            this.name = name; %
            this.stepper = stepper; %
            this.type = type; %

            this.existsType(type); % check if type is available

            this.content = cell(this.stepper.numsteps+1, 1);
        end
    end

    %% getdata %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        %> Returns the (coordinate) data on the grid of a step wrt the underlying
        %> basis (eg one vector for P0 or RT0, a matrix for P2P2, see Variable.type).
        %> @param this Instance of Variable.
        %> @param varargin arg1: Step number, data will be returned at @f$ t_\mathrm{stepno}@f$, arg2 (opt): function space
        %> @retval ret Numerical data set wrt the underlying basis. A vector-valued
        %>   variable may return a matrix with columns for each component.
        function ret = getdata(this, varargin)
            msg = 'HyPHM: Wrong syntax, please see documentation.';
            assert(nargin == 2 || nargin == 3, msg)
            stepno = varargin{1};

            % checking input arguments
            this.checkStep(this.stepper, stepno)
            if strcmp(this.name, ''), name = '[unnamed]';
            else name = this.name;
            end %#ok<PROP>
            assert(~isempty(this.content{stepno + 1, 1}), 'HyPHM: No data for Variable %s at step %d', name, stepno)

            % fetch data
            ret = this.content{stepno+1, 1};

            if nargin == 3 % transformation to different function space
                newtype = varargin{2};
                try
                    ret = eval(sprintf('%s.%sto%sslice(this.grid, ret)', this.type, this.type, newtype));
                catch
                    error('HyPHM: The conversion from %s to %s is not implemented, yet, or crashes.', this.type, newtype)
                end
            end
        end
    end

    %% distance %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        %> @brief Distance @f$\| w_h - w \|@f$ in a specified norm; note that by setting @f$w=0@f$ this method realizes norms @f$\| w_h\|@f$ .
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %> Evaluates the distance between the instance of
        %> <code>Variable</code> @f$w_h@f$ and a given algebraic function @f$w@f$ @f$\| w_h - w \|@f$
        %> at some specified time step, where
        %> @f$\|\cdot\|@f$ denotes the norm specified by one of the following tags:
        %>
        %> <table>
        %> <tr>
        %> <th align="left"><b>Norm</b></th>
        %> <th align="left"><b>Tag</b></th>
        %> <th align="left"><b>Description</b></th>
        %> </tr>
        %> <tr>
        %>  <td>@f$L^1(\Omega)@f$</td>
        %>  <td>@c L1</td>
        %>  <td>@f$\|v\|_{L^1(\Omega)} =   \int_\Omega |v| @f$ </td>
        %> </tr>
        %> <tr>
        %>  <td>@f$L^2(\Omega)@f$</td>
        %>  <td>@c L2</td>
        %>  <td>@f$\|v\|_{L^2(\Omega)} = \left( \int_\Omega |v|^2\right)^{1/2}@f$ </td>
        %> </tr>
        %> <tr>
        %>  <td>@f$\mathbb{P}_0(\mathcal{T})@f$</td>
        %>  <td>@c P0</td>
        %>  <td>@f$\|v\|_{\mathbb{P}_0(\mathcal{T})} = \left( \sum_{T\in\mathcal{T}}|T|\, v(x^\text{bary}_T)^2\right)^{1/2}@f$ (only for Variables of type @c P0) </td>
        %> </tr>
        %> <tr>
        %>  <td>@f$\mathbb{RT}_0(\mathcal{T})@f$</td>
        %>  <td>@c RT0</td>
        %>  <td>@f$\|\vec{v}\|_{\mathbb{RT}_0(\mathcal{T})} = \left( \sum_{T\in\mathcal{T}}\frac{|T|}{3}\sum_{E\subset T} \Big( \vec{v}(x^\text{bary}_E)\cdot\vec{\nu}_E\Big)^2\right)^{1/2}@f$ (only for Variables of type @c RT0) </td>
        %> </tr>
        %> </table>
        %>
        %> If @f$v@f$ is vector-valued, @f$|\cdot|@f$ denotes the Euclidean norm,
        %> ie @f$|\vec{v}| = \sqrt{v_1^2+v_2^2}@f$. The norm is obtained by
        %> numerical quadrature of high order.  This causes this
        %> method to be likely slow.
        %>
        %> @b Example:
        %> @code
        %> pres = Variable(g, st, 'pressure', 'P1');
        %> % [ computaion of pres ]
        %> presexact = @(t, x) 60*x(1)^2*x(2) - 20*x(2)^3;
        %> err = pres.distance(1, presexact, 'L2'); % the L2-discretization error
        %> @endcode
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %> @param this Instance of Variable.
        %> @param stepno Step number [1 int].
        %> @param fun  Algebraic function @f$w@f$ of the form <tt>f = @(t, x) f(t, x)</tt> [function_handle] where @f$\vec{x}\in\Omega@f$ [2x1].
        %> @param tag Tag which specifies the norm [string].
        %> @retval ret Distance of @f$w_h@f$ and @f$w@f$ in some norm, ie @f$\| w_h - w \|@f$.
        function ret = distance(this, stepno, fun, tag)
            assert(isa(fun, 'function_handle'), 'HyPHM: fun has to be a function_handle (t, x) -> IR');
            assert(ischar(tag), 'HyPHM: tag has to be a string, eg ''L1'' or ''L2''.');
            assert(strcmp(tag, 'P0') || strcmp(tag, 'RT0') || strcmp(tag, 'L1') || strcmp(tag, 'L2'), 'HyPHM: Tag has to be either L1 or L2')

            whdata = this.getdata(stepno); %#ok<NASGU>
            % fun is avaluated for a fixed t: statfun = statfun(x)
            statfun = @(x) fun(this.stepper.timeofstep(stepno), x); %#ok<NASGU>

            try
                ret = eval(sprintf('%s.distance%s(this.grid, whdata, statfun)', this.type, tag));
            catch
                error('HyPHM: The %s-distance computation for %s-functions is not implemented, yet, or crashes.', tag, this.type)
            end
        end
    end

    %% mean %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        %> @brief Integral mean @f$\frac{1}{|\Omega|}\int_\Omega w@f$.
        %>
        %> Evaluates the <i>integral mean</i>
        %> @f[\frac{1}{|\Omega|}\int_\Omega w(\vec{x})\, d\vec{x}@f]
        %> of a Variable @f$w@f$ by numerical quadrature for a given step number @c stepno.
        %> The integration is <i>exact</i> for the piecewise polynomials stored
        %> in Variable.
        %> For a vector-valued variable @f$\vec{w}@f$, the mean is computed
        %> <i>component-wise</i>, i.e.,
        %> @f$\mathtt{intmean}(\vec{w})\in\mathbb{R}^2@f$.
        %>
        %> @param this Instance of Variable.
        %> @param stepno Step number.
        %> @retval ret @f$\frac{1}{|\Omega|}\int_\Omega w@f$.
        function ret = mean(this, stepno)
            data = this.getdata(stepno); %#ok<NASGU>
            fprintf('%s.intmean(this.grid, data)', this.type)
            ret = eval(sprintf('%s.intmean(this.grid, data)', this.type));
        end
    end

    %% setdata %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        %> @brief Write data wrt the underlying discrete function space to
        %>        Variable.content.  This function examines the type of the
        %>        parameter @c data and passes it to the respective @c checkdata,
        %>        @c setfh2data or @c setfh3data
        %>        function in the package directory with the same name as
        %>        AbstractVariable.type.  If @c setfh2data or @c setfh3data is
        %>        called, the data is re-checked by ...
        %> @param this Instance of Variable.
        %> @param varargin (data) numerical data or function handle for all time points.
        %> @param varargin (stepno, data) ... only for time @f$ t_\mathrm{stepno}@f$.
        %> @todo  Postcheck data after setfun2data (over checkdata).
        %> @todo  Parallelize this!!!
        function setdata(this, varargin) % setdata(data) or setdata(stepno, data)
            msg = 'HyPHM: setdata(data) or setdata(stepno, data) where data is numerical or a function handle.';

            if nargin == 2 % do it for all time points
                if strcmp(this.name, '')
                    printline(3, 'Initializing an unnamed Variable for all time points...', this.name)
                else
                    %  printline(3, 'Initializing Variable ''%s'' for all time points...', this.name)
                end
                for k = 0:this.stepper.numsteps
                    this.setdata(k, varargin{:});
                end
            elseif nargin == 3 % set data for one single time point or a list of time points
                stepnums = varargin{1};
                data = varargin{2};
                for stepno = stepnums(:)' % loop over list of step numbers
                    this.checkStep(this.stepper, stepno)
                    if isnumeric(data)
                        eval(sprintf('%s.checkdata(this.grid, data)', this.type)); % verify data
                        this.content{stepno+1, 1} = data; % mount data
                    elseif isa(data, 'function_handle')
                        time = this.stepper.timeofstep(stepno); %#ok<NASGU>
                        if nargin(data) == 2
                            discdata = eval(sprintf('%s.setfh2data(time, this.grid, data)', this.type)); % compute data set via setfh2data
                            eval(sprintf('%s.checkdata(this.grid, discdata)', this.type)); % verify data
                            this.content{stepno+1, 1} = discdata; % mount data
                        elseif nargin(data) == 3
                            discdata = eval(sprintf('%s.setfh3data(time, this.grid, data)', this.type));
                            eval(sprintf('%s.checkdata(this.grid, discdata)', this.type)); % verify data
                            this.content{stepno+1, 1} = discdata; % mount data
                        else
                            error('HyPHM: Function handle has to be @(t, x) or @(T, X, Y).');
                        end
                    else
                        error(msg)
                    end
                end
            else % nargin ~= {2, 3}
                error(msg)
            end
        end

        %set data of Variable and save associated grid as well
        %data need to be given in vector format
        %to reduce memory consuption, instances of grid are reduced to most
        %frequently used properties
        function setdataGrid(this, stepno, data, grid)
            this.content{stepno+1, 1} = data;
            this.grid = grid;
            grid.baryE0T = 0;
            grid.V2T = 0;
            grid.V2E = 0;
            grid.A = 0;
            grid.b = 0;
            grid.idE = 0;
            grid.coordV0T = 0;
            grid.vecE = 0;
            grid.V0E = 0;
            grid.T0E = 0;
            grid.sigE0T = 0;
            grid.areaT = single(grid.areaT);
            grid.areaE = single(grid.areaE);
            grid.coordV = single(grid.coordV);
            grid.nuE = single(grid.nuE);
            this.grids{stepno+1, 1} = grid;
        end
    end

    %% visualize %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Generate a list of vtk-files (Doxygen documentation in visualize.m).
    methods % header
        visualize(this, varargin)
    end

    %% existsType %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %> @brief Check if a subpackage within +Variables with name 'type' exists (Static, Hidden).
    methods (Static, Hidden)
        function existsType(type)
            w = what('classes/Variables');
            z = w.packages;
            if ~any(ismember(z, type))
                printline(2, 'There are the following Variable types:')
                dir ./classes/Variables/+*
                error('HyPHM: Type ''%s'' not implemented.', type)
            end
        end
    end

    %% checkType %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %> @brief Check if a variable is of desired type.
    %> @b Example:
    %> @code
    %> variable.checkType('P0')
    %> variable.checkType({'P2P2', 'RT0'})
    %> @endcode
    %> @param this Instance of Variable.
    %> @param types String or cell of strings describing the valid discrete function spaces.
    methods (Hidden)
        function checkType(this, types)
            switch this.type
                case types
                    % do nothing
                otherwise
                    printline(-1, 'HyPHM: A Variable of the following type(s) is required:')
                    display(types)
                    error('HyPHM: See output above.')
            end
        end
    end

    %% checkStep %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Static, Hidden)
        %> Check whether the requested step number is in range.
        function checkStep(stepper, stepno)
            assert(stepno ~= -1, ...
                'HyPHM: Step number was -1, maybe you forgot to iterate the Stepper before solving: Stepper.next.')
            assert(stepno >= 0 && stepno <= stepper.numsteps, ...
                'HyPHM: Step number out of range [0, %d], requested step number was %d.', stepper.numsteps, stepno)
        end
    end

    %% SET %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        %> @brief Exclusively used by constructor.
        function set.grid(this, grid)
            msg = 'HyPHM: GRID has to be a instance of class Grid.';
            assert(isa(grid, 'AbstractGrid'), msg)
            this.grid = grid;
        end
        %> @brief Exclusively used by constructor.
        function set.stepper(this, stepper)
            msg = 'HyPHM: STEPPER has to be a instance of class Stepper.';
            assert(isa(stepper, 'Stepper'), msg)
            this.stepper = stepper;
        end
        %> @brief Exclusively used by constructor.
        function set.name(this, name)
            msg = 'HyPHM: NAME has to be a string.';
            assert(ischar(name), msg)
            this.name = name;
        end
        %> @brief Exclusively used by constructor.
        function set.type(this, type)
            msg = 'HyPHM: TYPE has to be a string.';
            assert(ischar(type), msg)
            this.type = type;
        end
    end

    methods (Static)
        function ret = isVariable(arg)
            ret = isa(arg, 'Variable');
        end
    end
end
