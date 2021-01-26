%> @file StokesLEVEL.m Class for the Stokes problem.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%>
%>   This class is an extension of Stokes, allowing for an additional level-set function argument L.
%>   All FEM in {L>0} will be disregarded, homogenous Dirichlet boundary conditions are
%>   automatically prescribed at {L=0}.
%>   An instance of this class should contain all data that are required by
%>   the solver Stokes.computeLevel(), i.e., coefficient functions, initial- and
%>   boundary data, contraints, and flags.  The coefficients and
%>   initial- and boundary data have to be given as instances of Variable (see documentation
%>   for the related discrete function spaces).
%>
%> @f[ \partial_t \vec{u} - \mu\,\nabla^2\vec{u} + \vec{\nabla} p = \vec{f}(t,\vec{x}) \qquad J\times\Omega@f]
%> @f[\vec{\nabla}\cdot\vec{u} = 0 \qquad J\times\Omega@f]
%> @f[ \vec{u} = \vec{u}^0   \qquad \left\{t_0\right\}\times\Omega@f]
%> @f[ \vec{u} = \vec{u}_\mathrm{D} \qquad J\times\Gamma_\mathrm{D}@f]
%> @f[ \big(\mu \mathbf{\nabla u} - p \mathbf{I}\big)\,\vec{\nu} = \vec{0} \qquad J\times\Gamma_\mathrm{N}@f]
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
%>     - StokesLEVEL.grid
%>     - StokesLEVEL.stepper
%>
%>
%>   <b>(ii) Unknowns</b>
%>
%>     - StokesLEVEL.U velocity @f$\vec{u}@f$ (type @c P2P2)
%>     - StokesLEVEL.P pressure @f$p@f$       (type @c P1)
%>
%>
%>   <b>(iii) Coefficient Functions</b>
%>
%>     - StokesLEVEL.N  dynamic viscosity @f$\mu@f$                 (type @c C)
%>     - StokesLEVEL.F  external force field  @f$\vec{f}@f$         (type @c P2P2)
%>     - StokesLEVEL.uD Dirichlet boundary @f$\vec{u}_\mathrm{D}@f$ (type @c P2P2, sampled at boundary only)
%>
%>
%>   <b>(iv) Flags</b>
%>
%>     - StokesLEVEL.isEvolution
%>
%>     This flag is set to false by default, ie the stationary Stokes problem is solved.
%>     To include the term @f$\partial_t \vec{u}@f$ set this flag to true.
%>
%>
%>   <b>(v) Constraints</b>
%>
%>     - StokesLEVEL.balanceP
%>
%>
%>   <b>(vi) Boundary Definition</b>
%>
%>     - StokesLEVEL.id2D
%>     - StokesLEVEL.id2N
%>     - StokesLEVEL.L    (implicit Dirichlet boundary)
%>
%> <h2>Methods</h2>
%>
%>   <b>(i) Constructor</b>
%>
%>     - StokesLEVEL(grid, stepper)
%>
%>   <b>(ii) Printing</b>
%>
%>     - StokesLEVEL.print()
%>
%>   <b>(iii) Solving</b>
%>
%>     - StokesLEVEL.computeLevel(Variable for @f$\vec{u}@f$, Variable for @f$p@f$)
%>
%>


classdef StokesLEVEL < AbstractProblem % < hgsetget & handle

    properties (SetAccess = private)
        %> Grid (set in construction)
        grid
        %> Stepper (set in construction)
        stepper
        %> name (set in construction)
        name
    end


    properties (Access = public)
        %> This may be defined if there's no Neumann bdry.  The balance only has effect on the
        %> pressure, not on the velocity solution.  If @c balanceP is equal to <tt>[]</tt> the
        %> contraint is not defined.
        balanceP

        %> If the equation contains the term @f$\partial_t \vec{u}@f$, default: false
        isEvolution
    end

    properties (Access = private, Hidden)
        globA
        globB
        globC
        globE
    end

    %% Unknowns %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % These are the unknowns which can be set from public.
    properties
        %> Velocity @f$\vec{u}@f$, Variable in @f$P_2(T)^2@f$.
        U
        %> Pressure @f$p@f$, Variable in @f$P_1(T)@f$.
        P
        %> Level-Set Function @f$L@f$, Variable in @f$P_1(T)@f$.
        L
    end

    %% Coefficients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % These are the coefficient functions which can be set from public.
    properties
        %> Dynamic viscosity @f$\nu@f$, Variable in @f$P_0(\Omega)@f$.
        N
        %> Source or sink @f$\vec{f}@f$, Variable in @f$P_2(T)^2@f$.
        F
        %> Dirichlet data @f$\vec{u}_\mathrm{D}@f$, Variable in @f$P_2(T)^2@f$.
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
        %> @param name Name of the problem.
        function this = StokesLEVEL(grid, stepper, name)
            msg = 'HyPHM: The constructor for this class (now) is StokesLEVEL(Grid, Stepper, name).';
            assert(nargin == 3, msg)
            assert(isa(grid, 'AbstractGrid'), msg)
            assert(isa(stepper, 'Stepper'), msg)
            assert(ischar(name))

            this.grid = grid;
            this.stepper = stepper;
            this.name = name;
            this.balanceP = []; % default: no contrains
            this.isEvolution = false; % default: no evolution term

            this.initVoidCoeffs('silent');
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

        function this = set.L(this, L)
            this.checkVariable(L)
            L.checkType('P1')
            this.L = L;
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

        function this = set.isEvolution(this, bool)
            assert(islogical(bool), 'HyPHM: StokesLEVEL.isEvolution expects true or false.')
            this.isEvolution = bool;
        end
    end

    %% Print Methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = public)
        %> Prints info acc to coefficients, flags and Newton.
        function print(this)
            this.printcoeffs
        end
    end

    % HEADERS
    %   NOTE: if you define headers, you cannot put the function file into the private directory!
    methods (Access = private, Hidden = true)
        printcoeffs(this)
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
