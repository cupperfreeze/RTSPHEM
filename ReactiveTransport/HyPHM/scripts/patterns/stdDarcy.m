%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @file stdDarcy.m
%> @brief HyPHM script for Darcy's law using the class Transport, see @ref stdDarcyTutorial.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @page stdDarcyTutorial Tutorial: The Darcy Problem
%>
%> @section s1 Problem Formulation
%>
%> Macroscopic saturated water flow in porous media.
%> Given a domain @f$\Omega@f$ with boundary
%> @f$\partial\Omega=\partial\Omega_\mathrm{D}\cup\partial\Omega_\mathrm{F}@f$ (see @ref notation),
%> we want to solve the stationary problem
%>
%> @f[\vec{\nabla}\cdot\vec{q} = f               \qquad \text{in}~\Omega @f]
%> @f[ \vec{q}=-\mathbf{K}\vec{\nabla}(\psi+z)   \qquad \text{in}~\Omega @f]
%> @f[ \psi = \psi_\mathrm{D}                    \qquad \text{on}~\partial\Omega_\mathrm{D}@f]
%> @f[ \vec{q}\cdot\vec{\nu} = q_\mathrm{F}      \qquad \text{on}~\partial\Omega_\mathrm{F}@f]
%>
%> with the unknowns @f$\vec{q}:\Omega\rightarrow\mathbb{R}^2@f$ and @f$\psi:\Omega\rightarrow\mathbb{R}@f$ denoting the <i>volumetric water flux density</i> and
%> the <i>pressure head</i>, respectively.  Furthermore, @f$\mathbf{K}:\Omega\rightarrow\mathbb{R}^{2,2}@f$
%> denotes the <i>hydraulic conductivity</i>, @f$z:\Omega\rightarrow\mathbb{R}@f$ the <i>elevation head</i> (the z-coordinate) and
%> @f$f:\Omega\rightarrow \mathbb{R}@f$ <i>sources</i> and/or <i>sinks</i>.  @f$\psi_\mathrm{D}:\partial\Omega_\mathrm{D}\rightarrow \mathbb{R}@f$ prescribes the pressure head on
%> the Dirichlet boundary, @f$g_\mathrm{flux}:\partial\Omega_\mathrm{flux}\rightarrow \mathbb{R}@f$ the flux through the boundary
%> (negative = inflow, positive = outflow, zero = no-flow).<br>
%>
%> This equation has the form of a diffusion/Poisson problem and can be
%> treated via the Transport class.<br>
%>
%> @section s2 The Domain
%> Let @f$\Omega=\,]0,\,2[\,\times\,]0,\,1[\,@f$ be our computational domain.
%> To generate rectangular bounded grids we can use the script
%> domainRectangle.m:
%> @code
%>   hmax = 2E-1; % mesh fineness
%>   g = domainRectangle(0, 2, 0, 1, hmax);
%>   g.print
%> @endcode
%> We want to check the edge IDs:
%> @code
%> g.visualize('idE')
%> @endcode
%> @image html stdDarcy-grid.png
%> Here, north edges are identified with ID 3, east with ID 2, etc and all
%> inner edges have always the ID zero.  This becomes important when the
%> boundary types have to be assigned.
%>
%> @section s3 The Time Adjustment
%> Although the problem is stationary, we need some time information.
%> We want to obtain the solution in one single step @f$t_0\rightarrow t_1@f$
%> where we arbitrarily choose @f$t_k= k@f$:
%> @code
%>   st = Stepper(0:1);
%>   st.print
%> @endcode
%>
%> @section s4 The Unknowns
%> The unknowns @f$\vec{q}@f$ and @f$\psi@f$ have to be defined via class
%> Variable.  Each instance of Variable lives in a discrete space which is
%> defined by the problem class, here Transport.  We have @f$\vec{q}\in\mathbb{RT}_0@f$
%> and @f$\psi\in \mathbb{P}_0@f$ (cf Variable.type):
%> @code
%>   phead = Variable(g, st, 'pressure head', 'P0');
%>   wflux = Variable(g, st, 'water flux', 'RT0');
%> @endcode
%> Further, the unknowns have to be initialized for @f$t_0@f$, ie step 0.
%> The way to do this is the method Variable.setdata.  One possibility is
%> @code
%>   phead.setdata(0, @(t, x) 0.0);
%>   wflux.setdata(0, @(t, x) [0.0; 0.0]);
%> @endcode
%>
%> @section s5 The Problem
%> With
%> @code
%>   darcy = Transport(g, st, 'Darcy problem');
%> @endcode
%> we obtain a problem with yet undefined coefficient functions and
%> boundary conditions.  Let's first define the <b>boundary types</b>:
%> Say we want the west and south boundary to be of inflow type, east of
%> no-flow type, ie impermeable boundary, north of type Dirichlet.  A
%> no-flow boundary is a flux boundary with @f$g_\mathrm{flux}=0@f$, so we
%> define
%> @code
%>   darcy.id2D = {3};
%>   darcy.id2F = {1,2,4};
%> @endcode
%> The <b>boundary data</b> has to be of class Variable, similar to the unknowns.
%> The required type can again be found in Transport.uD and Transport.gF.
%> Lets choose
%> @code
%>   darcy.uD.setdata(@(t, x) 1.0);
%> @endcode
%> Analogeous,
%> @code
%>   darcy.gF.setdata(@(t, x) x(1)-2);
%> @endcode
%> We now have to define the <b>coefficient functions</b> such that the
%> transport equation, simplifies to the Darcy problem.
%> Comparing both equations, we see that the coefficients A, B, C, E in
%> Transport have to vanish.  We do so by zero-initialization or not
%> defining them.  Remaining the diffusion coefficient and the source term.
%> Comparing the definition of the boundary data above, we can alternatively
%> mount the instances of Variable directly into the problem:
%> @code
%> darcy.D.setdata(eye(2)); % D is a 2 times 2 matrix, all 4 entries should be constant (C)
%> darcy.F.setdata(@(t, x) 2*(1-norm(x-[1;0.5])));
%> @endcode
%> If we are unsure and want to check the definitions, we can visualize
%> them by
%> @code
%> darcy.F.visualize
%> @endcode
%> @image html stdDarcySource.jpg
%>
%> @section s6 Assembly and Computation
%> The Stepper @c st points to step zero at the moment.  To obtain the
%> solution of the following step, step one, we have to iterate it.
%> Thereafter we call the solver:
%> @code
%> st.next;
%> darcy.computeLevel;
%> @endcode
%>
%> @section s7 Visualization
%> Every instance of Variable can be visualized if it contains data for
%> every time step.  The raw data can be obtained via Variable.getdata.
%> @code
%> phead.visualize
%> wflux.visualize
%> @endcode
%> The @ref paraview output may look like
%> @image html stdDarcySolution.jpg
%>
%> @section complete The Complete Script
%> All together, the script may read
%> @code
%> %% Grid
%> hmax = 2E-1; % mesh fineness
%> g = domainRectangle(0,2, 0,1, hmax);
%> g.print
%> g.visualize('idE')
%>
%> %% Stepper
%> st = Stepper(0:1);
%> st.print
%>
%> %% Unknowns
%> phead = Variable(g, st, 'pressure head', 'P0');
%> wflux = Variable(g, st, 'water flux',    'RT0');
%> phead.setdata(0, @(t, x) 0.0);
%> wflux.setdata(0, @(t, x) [0.0; 0.0]);
%>
%> %% Problem Data
%> darcy = Transport(g, st, 'Darcy problem');
%> % unknowns
%> darcy.Q = wflux;
%> darcy.U = phead;
%> % flags
%> darcy.isDstationary = true;
%> % boundary data
%> darcy.id2D = {3};
%> darcy.id2F = {1,2,4};
%> darcy.uD.setdata(@(t, x) 1.0);
%> darcy.gF.setdata(@(t, x) x(1)-2);
%> % coefficient functions
%> darcy.D.setdata(eye(2));
%> darcy.F.setdata(@(t, x) 2*(1-norm(x-[1;0.5])));
%>
%> %% Assembly and Computation
%> st.next;
%> darcy.computeLevel;
%>

%% Visualization
%> darcy.F.visualize % source/sink term
%> wflux.visualize   % water flux
%> phead.visualize   % pressure head
%> @endcode
%>
%>
%> @section s8 Possible Extensions
%> Rather to use a constant coefficient
%>@f$\mathbf{K}:\Omega\rightarrow\mathbb{R}^{2,2}@f$,
%>you may want to have a spatially varying one.  Since <code>darcy.D</code>
%>is of type <code>CCCC</code> (constant in each component) by default (cf. <code>Transport.D</code>),
%>we first have to create a @c Variable of type <code>P0P0P0P0</code> (elementwise
%>constant in each component, see <code>Variable.type</code>):
%>@code
%>conduct = Variable(g, st, 'conductivity', 'P0P0P0P0');
%>darcy.D = conduct;
%>@endcode
%>or simply
%>@code
%>darcy.D = Variable(g, st, 'conductivity', 'P0P0P0P0');
%>@endcode
%>This variable may be initialized using the <code>@(t, x)</code>-syntax, where the time @f$t@f$ is irrelevant here.  For instance,
%>if we want to define
%>@f[ \mathbf{K}(\vec{x}) = \begin{cases}0.1&\text{if}~x\in[0.25,0.75]\\1&\text{otherwise}\end{cases} @f]
%>we write
%>@code
%>darcy.D.setdata(@(t, x) eye(2)*(1 - 0.9*(x(1)<=0.75 && x(1)>=0.25)));
%>@endcode
%>The expression that defines the <code>function_handle</code> does not contain the @c t in the argument.
%>However, this is the syntax expected from @c Variable.setdata.
%>
%>The result looks like
%> @image html stdDarcy-sol2.png
%> where we used 'Stream Tracer' in Paraview to visualize the variable
%> <code>wflux</code> after calling the filter 'Delaunay 2D'.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Grid
hmax = 2E-1; % mesh fineness
g = domainRectangle(0, 2, 0, 1, hmax);
g.print
% g.visualize('idE')

%% Stepper
st = Stepper(0:1);
st.print

%% Unknowns
phead = Variable(g, st, 'pressure head', 'P0');
wflux = Variable(g, st, 'water flux', 'RT0');
phead.setdata(0, @(t, x) 0.0);
wflux.setdata(0, @(t, x) [0.0; 0.0]);

%% Problem Data
darcy = Transport(g, st, 'Darcy problem');
% unknowns
darcy.Q = wflux;
darcy.U = phead;
% boundary data
darcy.id2D = {3};
darcy.id2F = {1, 2, 4};
darcy.uD.setdata(@(t, x) 1.0);
darcy.gF.setdata(@(t, x) x(1)-2);
% coefficient functions
darcy.F.setdata(@(t, x) 1*(1 - norm(x - [0.1; 0.5])));

%% Assembly and Computation
st.next;
darcy.computeLevel;

%% Visualization
darcy.F.visualize % source/sink term
wflux.visualize % water flux
phead.visualize % pressure head
