%> @file Perm.m Class for permeability cell problems.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%>
%>   An instance of this class is only supposed to have a level-set function L initialized.  It  is to be given as instances of Variable (see documentation
%>   for the related discrete function spaces). This is an optimized routine for efficient computation of permeability tensors. It is built upon the class Stokes.m but it is tailored to that specific use.
%>   All coefficients are automatically initialized with their correct value. Accordingly, the interior boundary is set to be a homogeneous Dirichlet boundary.
%>
%> @f[ \nabla^2\vec{u_i} + \vec{\nabla} p = \vec{e}_i \qquad \Omega@f]
%> @f[\vec{\nabla}\cdot\vec{u_i} = 0 \qquad \Omega@f]
%> @f[ \vec{u} = 0 \qquad \Gamma_\mathrm{int}@f]
%> @f[i =1,2 @f]
%>
%>
%> <h2>Properties</h2>
%>
%>   <b>(i) Grid and Stepper</b>
%>
%>     - Perm.grid
%>     - Perm.stepper
%>
%>
%>   <b>(ii) Unknowns</b>
%>
%>     - Perm.U  velocity @f$\vec{u}@f$  (type @c P2P2 or P1BubbleP1Bubble)
%>     - Perm.U2 velocity @f$\vec{u2}@f$ (type @c P2P2 or P1BubbleP1Bubble)
%>     - Perm.P  pressure @f$p@f$        (type @c P1)
%>
%>
%>   <b>(iii) Coefficient Functions</b>
%>
%>     - Perm.L  Level-Set Function @f$L@f$                 (type @c P1)
%>
%>
%>   <b>(iv) Flags</b>
%>
%>     - none needed
%>
%>   <b>(vi) Boundary Definition</b>
%>
%>     - none needed
%>
%>
%> <h2>Methods</h2>
%>
%>   <b>(i) Constructor</b>
%>
%>     - Perm(grid, stepper)
%>
%>   <b>(ii) Printing</b>
%>
%>     - Perm.print()
%>
%>   <b>(iii) Solving</b>
%>
%>     - Perm.computeLevel(Variable for @f$\vec{u}@f$, Variable for @f$p@f$)
%>
%>

classdef Perm < AbstractProblem % < hgsetget & handle

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
        %> Velocity @f$\vec{u}@f$, Variable in @f$P_2(T)^2@f or @f$P_1Bubble(T)^2@f$.
        U
        U2
        %> Pressure @f$p@f$, Variable in @f$P_1(T)@f$.
        P
        L %LevelSet
    end

    %% Coefficients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % These are the coefficient functions which can be set from public.
    properties
        %> Dynamic viscosity @f$\nu@f$, Variable in @f$P_0(\Omega)@f$.
        N
        %> Source or sink @f$\vec{f}@f$, Variable in @f$P_2(T)^2@f$.
        %F
        %> Dirichlet data @f$\vec{u}_\mathrm{D}@f$, Variable in @f$P_2(T)^2@f$.
        uD
        %> True if P1Bubble/P1 discretization in use
        Bubble
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
        function this = Perm(grid, stepper, name, varargin)
            msg = 'HyPHM: The constructor for this class (now) is Perm(Grid, Stepper, name).';
            assert(nargin >= 3, msg)
            assert(isa(grid, 'AbstractGrid'), msg)
            assert(isa(stepper, 'Stepper'), msg)
            assert(ischar(name))

            this.grid = grid;
            this.stepper = stepper;
            this.name = name;
            this.balanceP = []; % default: no contrains
            this.isEvolution = false; % default: no evolution term

            Bubble = false;
            if ismember(varargin, 'Bubble')
                Bubble = true;
            end

            this.Bubble = Bubble;
            this.initVoidCoeffs('silent');
        end % constructor
    end

    %% Set Methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        function this = set.U(this, U)
            this.U = U;
        end

        function this = set.U2(this, U2)
            this.U2 = U2;
        end

        function this = set.P(this, P)
            this.P = P;
        end

        function this = set.L(this, L)
            this.checkVariable(L)
            L.checkType('P1')
            this.L = L;
        end

        function this = set.N(this, N)
            this.N = N;
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
