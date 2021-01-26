%> @file stdNavierStokes.m Standard script for solving the Navier-Stokes problem.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> Dependencies:
%>   - navStokes_F.m
%>   - navStokes_uD.m


addpath('./scripts/scripts-patterns/stdNavierStokes') % add folder to path where scenario stuff is located

%% Initialization of Timer
st = Stepper(0:.1:.5);
st.status % status bar
st.print % print stepper data to command window

%% Initialization of Grids %%
hmax = 2E-1;
g = domainRectangle(-1, 1, -1, 1, hmax);
% g.visualize('idE')
g.print

%% Solution Vectors
pres = Variable(g, st, 'pressure', 'P1');
velo = Variable(g, st, 'velocity', 'P2P2');

%% Initialization of Problem Data
d = NavierStokes(g, st);

% Flags
d.isConvection = true; % convective term?
d.isEvolution = true; % evolution term?
d.isSymmetric = true; % which kind of diffusion tensor?

% Constraints
d.balanceP = []; % (default)    % mean pressure equal to a value? (only if there's no Neumann bdry)

% Parameters and Coefficient Functions
d.N = Variable(g, st, 'viscosity', 'C');
d.N.setdata(1E-2);
d.F = Variable(g, st, 'source', 'P2P2');
d.F.setdata(@navStokes_F);
d.uD = Variable(g, st, 'bdrydata', 'P2P2');
d.uD.setdata(@navStokes_uD);

% Boundary IDs
d.id2D = {2, 3, 4};
d.id2N = {1};

% Parameters for the Newton iteration
% d.maxIter = 20;
% d.maxRes  = 1E-10;

%% Initialization of Initial Solution
% initializaton not needed for discretization, though for visualization
pres.setdata(0, zeros(g.numV, 1));
velo.setdata(0, zeros(g.numV + g.numE, 2));

%% Assembly and Computation
printline(1, 'Assembly and Computation');
while st.next
    d.computeLevel(velo, pres);
end % end stepper iteration

%% Visualization
printline(1, 'Visualization');
pres.visualize
velo.visualize
