%> @file NavierStokes.m Class for the <b>Navier&ndash;Stokes</b> or the <b>instationary Stokes</b> problem.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @brief Class for the <b>Navier&ndash;Stokes</b> or the <b>instationary Stokes</b> problem.
%>
%>   An instance of this class should contain all data,  which are required by
%>   solver(), i.e., coefficient functions, initial- and boundary data,
%>   reaction parameters, species properties an so on.  The coefficients and
%>   initial- and boundary data have to be given as instances of Variable.
%>
%> @f[ \partial_t (A \vec{u}) - \vec{\nabla}\cdot\big(\nu (\vec{\nabla}\vec{u} + {\vec{\nabla}\vec{u}}^\mathrm{T})\big) + \left(\vec{u}\cdot\vec{\nabla}\right)\vec{u} + \vec{\nabla} p = \vec{F} \qquad J\times\Omega@f]
%> @f[\vec{\nabla}\cdot\vec{u} = 0 \qquad J\times\Omega@f]
%> @f[ \vec{u} = \vec{u}^0   \qquad \left\{t_0\right\}\times\Omega@f]
%> @f[ \vec{u} = \vec{u}_\mathrm{D} \qquad J\times\Gamma_\mathrm{D}@f]
%> @f[ \big(\nu (\vec{\nabla}\vec{u} + {\vec{\nabla}\vec{u}}^\mathrm{T}) - p \mathbf{I}\big)\cdot\vec{\nu} = \vec{0} \qquad J\times\Gamma_\mathrm{N}@f]
%>
%>
%>
%>
%> In general, the coefficient functions have to be mounted and the flags to be modified
%> according to the desired problem.  Moreover, one has to specify the sets
%> of edge ids for the desired boundary types.
%>
%>
%> <h2>Properties</h2>
%>
%>   <b>(i) Grid and Stepper</b>
%>
%>     - NavierStokes.grid
%>     - NavierStokes.stepper
%>
%>
%>   <b>(ii) Unknowns</b>
%>
%>     - NavierStokes.U velocity @f$\vec{u}@f$ (type @c P2P2)
%>     - NavierStokes.P pressure @f$p@f$       (type @c P1)
%>
%>
%>  <b>(iii) Coefficient Functions</b>
%>
%>     - NavierStokes.N
%>     - NavierStokes.F
%>     - NavierStokes.uD
%>   <b>(iv) Flags</b>
%>
%>     - NavierStokes.isConvection
%>     - NavierStokes.isEvolution
%>     - NavierStokes.isSymmetric
%>
%>     All flags are set to false by default.
%>
%>
%>   <b>(v) Constraints</b>
%>
%>     - NavierStokes.balanceP
%>
%>
%>   <b>(vi) Boundary Definition</b>
%>
%>     - NavierStokes.id2D
%>     - NavierStokes.id2N
%>
%>
%>
%>   <b>(vii) Parameter for the Newton Solver (convection term)</b>
%>     maxIter
%>     maxRes
%>
%>
%> <h2>Methods</h2>
%>
%>   <b>(i) Constructor</b>
%>
%>     - NavierStokes(grid, stepper)
%>
%>   <b>(ii) Printing</b>
%>
%>     - NavierStokes.print()
%>
%>   <b>(iii) Solving</b>
%>
%>     - NavierStokes.computeLevel(Variable for @f$\vec{u}@f$, Variable for @f$p@f$)
%>
%>


classdef NavierStokes < AbstractProblem % < hgsetget & handle

    properties (SetAccess = private)
        %> Grid (set in construction)
        grid
        %> Stepper (set in construction)
        stepper
    end


    properties (Access = public)
        %> @brief Flag for non-linear convection term.
        %> If unset, this is the Stokes equation [logical].
        isConvection
        %> @brief Flag for temporal derivative term [logical].
        %> If unset, this is the stationary version.
        isEvolution
        %> @brief Flag for the form of the diffusion tensor [logical].
        %> If set true, diffusion is @f$\nabla\cdot(\nu(\nabla u + (\nabla u)^T) )@f$,
        %> otherwise @f$\nabla\cdot(\nu(\nabla u))@f$ (default).
        isSymmetric
        %> @brief Pressure balance contraint, ie mean pressure equals @c balanceP [logical].
        %>
        %> This may be defined if there's no Neumann bdry.  The balance only has effect on the
        %> pressure, not on the velocity solution.  If @c balanceP is equal to <tt>[]</tt> the
        %> contraint is not defined.
        balanceP
    end

    % Parameters for Newton's method (set to default in constructor)
    properties (Access = public)
        %> Maximum Newton iterates  [scalar].
        maxIter
        %> Residual bound for Newton [scalar].
        maxRes
    end

    %% Unknowns %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % These are the unknowns which can be set from public.
    properties
        %> Velocity @f$\vec{u}@f$, Variable in @f$P_2(T)^2@f$.
        U
        %> Pressure @f$p@f$, Variable in @f$P_1(T)@f$.
        P
    end

    %% Coefficients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % These are the coefficient functions which can be set from public.
    properties
        %> Dynamic viscosity @f$\nu@f$, Variable in @f$P_0(\Omega)@f$.
        N
        %> Source or sink, Variable in @f$P_2(T)^2@f$.
        F
        %> Dirichlet data, Variable in @f$P_2(T)^2@f$.
        uD
    end

    % Edge type definition
    properties (Access = public)
        %> @brief Set of edge ids which should identify Neumann edges [cell].
        id2N
        %> @brief Set of edge ids which should identify Dirichlet edges [cell].
        %> Example:
        %> If, eg, edges with ids 2 and 5 are desired to become Dirichlet, set <tt>id2D = {2, 5}</tt>.
        id2D
    end

    %% Constructor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        %> @brief Constructor.
        %> @param grid Instance of Grid.
        %> @param stepper Instance of Stepper.
        function this = NavierStokes(grid, stepper)
            msg = 'HyPHM: Constructor is NavierStokes(Grid, Stepper).';
            assert(nargin == 2, msg)
            assert(isa(grid, 'AbstractGrid'), msg)
            this.grid = grid;
            assert(isa(stepper, 'Stepper'), msg)
            this.stepper = stepper;

            this.isConvection = false;
            this.isEvolution = false;
            this.isSymmetric = false;

            this.balanceP = [];

            this.maxRes = 1E-3;
            this.maxIter = 20;
        end % constructor
    end

    %% Set Methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        function this = set.U(this, U)
            this.checkVariable(U)
            U.checkType('P2P2')
            this.U = U;
        end

        function this = set.P(this, P)
            this.checkVariable(P)
            P.checkType('P1')
            this.P = P;
        end

        function this = set.N(this, N)
            this.checkVariable(N)
            N.checkType('C')
            this.N = N;
        end

        function this = set.F(this, F)
            this.checkVariable(F)
            F.checkType('P2P2')
            this.F = F;
        end

        function this = set.uD(this, uD)
            this.checkVariable(uD)
            uD.checkType('P2P2')
            this.uD = uD;
        end
    end

    %% Print Methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = public)
        %> Prints info acc to coefficients, flags and Newton.
        function print(this)
            this.printnewton
            this.printflags
        end
    end

    % HEADERS
    %   NOTE: if you define headers, you cannot put the function file into the private directory!
    methods (Access = private, Hidden = true)
        printnewton(this)
        printflags(this)
    end

    %% Solvers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        solve(this, U, P, varargin)
    end

    methods (Static, Hidden)

        %% checkVariable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %> @brief Check if argument is an instance of Variable.
        function checkVariable(arg)
            assert(isa(arg, 'Variable'), ...
                'HyPHM: Instance of Variable required.')
        end
    end
end
