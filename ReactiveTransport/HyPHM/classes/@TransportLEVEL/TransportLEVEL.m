%> @file TransportLEVEL.m Class for the transport problem.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> This is a specialized routine to solve transport problems on a (time-dependent) domain
%> characterized by a level-set function L. The resulting interior boundary is automatically
%> identified as a flux boundary.
%> Any triangle with all vertices Level-set value >=-eps are regarded as
%> solid and do not belong to fluid domain.
%> An instance of this class should contain all data that are required by
%> TransportLEVEL.computeLevel(), i.e., coefficient functions, initial, and boundary data,
%> constraints and so on.  The coefficients and the data have to be given as instances of Variable.
%>
%> @section sec1 Problem
%> for the following, we set \Omega = \Omega \cap {L<0}
%> @f[ \vec{q} = -\mathbf{D}\vec{\nabla}u + \vec{C}u + \vec{E}  \qquad J\times\Omega @f]
%> @f[ \partial_t (Au + B) + \vec{\nabla}\cdot \vec{q} = F  \qquad J\times\Omega @f]
%> @f[ u = u^0   \qquad \left\{t_0\right\}\times\Omega@f]
%> @f[ u = u_\mathrm{D} \qquad J\times\partial\Omega_\mathrm{D}@f]
%> @f[ (-\mathbf{D}\vec{\nabla} u)\cdot\vec{\nu} = 0 \qquad J\times\partial\Omega_\mathrm{N}@f]
%> @f[ \vec{q}\cdot\vec{\nu} = g_\mathrm{flux} \qquad J\times\partial\Omega_\mathrm{flux}@f]
%>
%> In general, the coefficient functions have to be mounted and the flags to be modified
%> according to the desired problem.  Moreover, one has to specify the sets
%> of edge IDs for the desired boundary types (see @ref stdDarcyTutorial).
%>
%> @section sec2 Unknowns and Coefficients
%> The coefficients will be initialized by default when a variable of
%> TransportLEVEL is generated (see also TransportLEVEL.initVoidCoeffs).  You can
%> either define new instances of Variable and mount those to the problem
%> (i.e. to an instance of TransportLEVEL) or modify them directly:
%>
%> @code
%> % st = Stepper(...);
%> % g  = Grid(...);
%> t = TransportLEVEL(s, st, 'my TransportLEVEL problem');
%> convection = Variable(g, st, 'convection', 'RT0');
%> convection.setdata(@(t, x) t*[x(1), -x(2)]);
%> t.C = convection;
%> @endcode
%>
%> or
%>
%> @code
%> % st = Stepper(...);
%> % g  = Grid(...);
%> t = TransportLEVEL(s, st, 'my TransportLEVEL problem');
%> t.C.setdata(@(t, x) t*[x(1), -x(2)]);
%> @endcode
%>
%> In the latter example, TransportLEVEL.C.name is still <code>'C'</code> while
%> in the first one, TransportLEVEL.C.name is <code>'convection'</code>.
%>
%> The default initial values for the unknowns and coefficients are:
%>
%> <table>
%> <tr>
%> <th align="center" width="100"><b>Quantity</b></th>
%> <th align="center" width="100"><b>Variable</b></th>
%> <th align="center" width="100"><b>Default Name</b></th>
%> <th align="center" width="100"><b>Space</b></th>
%> <th align="center" width="250"><b>Default Value</b></th>
%> </tr>
%> <tr>
%>  <td>@f$\vec{q}@f$</td>
%>  <td>TransportLEVEL.Q</td>
%>  <td>@c Q</td>
%>  <td>@f$\mathbb{RT}_0(\mathcal{T})@f$</td>
%>  <td>@f$0@f$ (only at @f${t_0}@f$, all other time slices remain uninitialized)</td>
%> </tr>
%> <tr>
%>  <td>@f$u@f$</td>
%>  <td>TransportLEVEL.U</td>
%>  <td>@c U</td>
%>  <td>@f$\mathbb{P}_0(\mathcal{T})@f$</td>
%>  <td>@f$0@f$ (only at @f${t_0}@f$, all other time slices remain uninitialized)</td>
%> </tr>
%> <tr>
%>  <td>@f$A@f$</td>
%>  <td>TransportLEVEL.A</td>
%>  <td>@c A</td>
%>  <td>@f$\mathbb{P}_0(\mathcal{T})@f$</td>
%>  <td>@f$0@f$</td>
%> </tr>
%> <tr>
%>  <td>@f$B@f$</td>
%>  <td>TransportLEVEL.B</td>
%>  <td>@c B</td>
%>  <td>@f$\mathbb{P}_0(\mathcal{T})@f$</td>
%>  <td>@f$0@f$</td>
%> </tr>
%> <tr>
%>  <td>@f$\vec{C}@f$</td>
%>  <td>TransportLEVEL.C</td>
%>  <td>@c C</td>
%>  <td>@f$\mathbb{RT}_0(\mathcal{T})@f$</td>
%>  <td>@f$\vec{0}@f$</td>
%> </tr>
%> <tr>
%>  <td>@f$\mathbf{D}@f$</td>
%>  <td>TransportLEVEL.D</td>
%>  <td>@c D</td>
%>  <td>@f$\mathbb{P}_0(\Omega)^{2, 2}@f$ (default) or @f$\mathbb{P}_0(\mathcal{T})^{2, 2}@f$</td>
%>  <td>@f$\mathbf{1}@f$ (unit matrix)</td>
%> </tr>
%> <tr>
%>  <td>@f$\vec{E}@f$</td>
%>  <td>TransportLEVEL.E</td>
%>  <td>@c E</td>
%>  <td>@f$\mathbb{RT}_0(\mathcal{T})@f$</td>
%>  <td>@f$\vec{0}@f$</td>
%> </tr>
%> <tr>
%>  <td>@f$F@f$</td>
%>  <td>TransportLEVEL.F</td>
%>  <td>@c F</td>
%>  <td>@f$\mathbb{P}_0(\mathcal{T})@f$</td>
%>  <td>@f$0@f$</td>
%> </tr>
%> <tr>
%>  <td>@f$g_\mathrm{flux}@f$</td>
%>  <td>TransportLEVEL.gF</td>
%>  <td>@c gF</td>
%>  <td>@f$\mathbb{P}_0(\mathcal{E})@f$</td>
%>  <td>@c NaN</td>
%> </tr>
%> <tr>
%>  <td>@f$u_\mathrm{D}@f$</td>
%>  <td>TransportLEVEL.uD</td>
%>  <td>@c uD</td>
%>  <td>@f$\mathbb{P}_0(\mathcal{E})@f$</td>
%>  <td>@c NaN</td>
%> </tr>
%> </table>
%>
%> The function spaces are listed in Variable.type.
%>
%> @section Upwinding
%> If the <i>local P&eacute;clet number</i>
%> @f[\mathrm{Pe}_T = \frac{\|\vec{C}\|_{L^\infty}(T)\,h_T}{2\|\mathbf{D}\|_{L^\infty}(T)}@f]
%> for triangles @f$T\in\mathcal{T}_h@f$ is large, we might need to stabilize the scheme in order to suppress
%> spurious oscillations.  We can do so by setting the flag
%> <code>TransportLEVEL.isUpwind</code>, which is <code>'none'</code> by default:

%> <table>
%> <tr>
%> <th align="center"><b>Upwinding Strategy</b></th>
%> <th align="center"><b>Tag</b></th>
%> <th align="center"><b>Comment</b></th>
%> </tr>
%> <tr>
%>  <td>no upwinding <i>(default)</i></td>
%>  <td><code>'none'</code></td>
%>  <td></td>
%> </tr>
%> <tr>
%>  <td>full upwinding</td>
%>  <td><code>'full'</code></td>
%>  <td></td>
%> </tr>
%> <tr>
%>  <td>exponential upwinding</td>
%>  <td><code>'exp'</code></td>
%>  <td bgcolor="#FFFFFF">@image html isUpwindExp.png</td>
%> </tr>
%> <tr>
%>  <td>exponential upwinding with switch</td>
%>  <td><code>'alt'</code></td>
%>  <td bgcolor="#FFFFFF">@image html isUpwindAlt.png</td>
%> </tr>
%> </table>

classdef TransportLEVEL < AbstractProblem

    properties (SetAccess = private)
        %> <i>(SetAccess = private)</i> Instance of Grid or FoldedGrid (set in construction).
        grid
        %> <i>(SetAccess = private)</i> Instance of Stepper (set in construction).
        stepper
        %> <i>(SetAccess = private)</i> Name of the problem [string] (set in construction).
        name
    end

    properties (Access = public)
        %> Balance constraints, see NavierStokes.
        balanceU
        %> Upwind flag: 'none' (<i>default</i>), 'full', 'exp' (exponential upwinding), 'alt' (exponential upwinding if local Pe number > 2)
        isUpwind
    end

    % Edge type definition (set to default in constructor)
    properties (Access = public)
        %> List of edge IDs defining Neumann boundaries [cell].
        id2N
        %> List of edge IDs defining flux boundaries [cell].
        id2F
        %> List of edge IDs defining Dirichlet boundaries [cell].
        id2D
    end

    % Storage of stationary assembling matrices
    properties (Access = private, Hidden)
        %> <i>(Access = private, Hidden)</i> Storage of global matrix @f$\mathbf{B}@f$, see @ref Frank2013, p.74.
        globB
        %> <i>(Access = private, Hidden)</i> Storage of global matrix @f$\mathbf{D}@f$, see @ref Frank2013, p.74.
        globD
        %> <i>(Access = private, Hidden)</i> Storage of global matrix @f$\mathbf{H}@f$, see @ref Frank2013, p.74.
        globH
    end

    %% Unknowns %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % These are the unknowns  which can be set from public.
    properties
        %> Mass flux in transport or water flux in Darcy, @f$\vec{q}@f$, Variable in @f$\mathbb{RT}_0(T)@f$.
        Q
        %> Concentration in transport or pressure head in Darcy, @f$u@f$, Variable in @f$\mathbb{P}_0(\mathcal{T})@f$.
        U
    end

    %% Coefficients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % These are the coefficient functions which can be set from public.
    properties
        %> Porosity in transport equation, Variable in @f$\mathbb{P}_0(\mathcal{T})@f$.
        A
        %> Porosity in pressure equation, Variable in @f$\mathbb{P}_0(\mathcal{T})@f$.
        B
        %> Convection, Variable in @f$\mathbb{RT}_0(\mathcal{T})@f$.
        C
        %> Diffusion tensor, Variable in @f$\mathbb{P}_0(\mathcal{T})^{2, 2}@f$ or @f$\mathbb{P}_0(\Omega)^{2, 2}@f$.
        D
        %> Force/Gravitation, Variable in @f$\mathbb{RT}_0(\mathcal{T})@f$.
        E
        %> Source or sink term, Variable in @f$\mathbb{P}_0(\mathcal{T})@f$.
        F
        %> Dirichlet data, Variable in @f$\mathbb{P}_0(\mathcal{E})@f$.
        uD
        %> Flux data, Variable in @f$\mathbb{P}_0(\mathcal{E})@f$.
        %> @todo allow C(E).
        gF
        %> Level-Set Function, Variable in @f$\mathbb{P}_1(\mathcal{T})@f$
        L
    end

    %% Constructor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        %> @brief Constructor.
        %> The balance contraints are set to empty, ie deactivated.
        %> @param grid Instance of Grid.
        %> @param stepper Instance of Stepper.
        %> @param name Name of the problem.
        function this = TransportLEVEL(grid, stepper, name)
            msg = 'HyPHM: The constructor for this class (now) is TransportLEVEL(Grid, Stepper, name).';
            assert(nargin == 3, msg)
            assert(isa(grid, 'AbstractGrid'), msg)
            assert(isa(stepper, 'Stepper'), msg)
            assert(ischar(name))

            this.grid = grid;
            this.stepper = stepper;
            this.name = name;
            this.balanceU = []; % no balance constraint at initialization
            this.isUpwind = 'none';

            this.initVoidCoeffs('silent');
        end % constructor
    end

    %% Set Methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods


        function this = set.Q(this, Q)
            this.checkVariable(Q)
            Q.checkType('RT0')
            this.Q = Q;
        end

        function this = set.U(this, U)
            this.checkVariable(U)
            U.checkType('P0')
            this.U = U;
        end

        function this = set.A(this, A)
            this.checkVariable(A)
            A.checkType('P0')
            this.A = A;
        end

        function this = set.B(this, B)
            this.checkVariable(B)
            B.checkType('P0')
            this.B = B;
        end

        function this = set.C(this, C)
            this.checkVariable(C)
            C.checkType('RT0')
            this.C = C;
        end

        function this = set.D(this, D)
            this.checkVariable(D)
            D.checkType({'P0P0P0P0', 'CCCC'})
            this.D = D;
        end

        function this = set.E(this, E)
            this.checkVariable(E)
            E.checkType('RT0')
            this.E = E;
        end

        function this = set.F(this, F)
            this.checkVariable(F)
            F.checkType('P0')
            this.F = F;
        end

        function this = set.uD(this, uD)
            this.checkVariable(uD)
            uD.checkType('P0E')
            this.uD = uD;
        end

        function this = set.L(this, L)
            this.checkVariable(L)
            L.checkType('P1')
            this.L = L;
        end

        function this = set.gF(this, gF)
            this.checkVariable(gF)
            gF.checkType('P0E')
            this.gF = gF;
        end
    end % set functions

    methods
        function print(this)
            % @todo Implement this!
        end
    end

    methods (Static, Hidden)

        %% checkVariable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %> @brief Check if argument is an instance of Variable.
        function checkVariable(arg)
            assert(isVariable(arg), 'HyPHM: Instance of Variable required.')
        end
    end


end
