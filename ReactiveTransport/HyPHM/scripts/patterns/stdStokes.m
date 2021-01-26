%> @file stdStokes.m Standard script for solving the Stokes problem.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization of Timer %%
st = Stepper(0:1);
st.status % status bar
st.print % print stepper data to command window

%% Initialization of Grids %%
% hmax = 5E-2;
hmax = 0.1;
g = domainRectangle(-1, 1, -1, 1, hmax);
g.print

%% Initialization of Problem Data
s = Stokes(g, st, 'Stokes');

% constraints
s.balanceP = []; % (default)    % mean pressure equal to a value? (only if there's no Neumann bdry)

% coefficients
s.uD.setdata(@(t, x) exp(t)*[cos(x(1)); sin(x(2))]);
s.F.setdata(@(t, x) [x(1); x(2)]);

% boundary IDs
s.id2D = {2, 3, 4};
s.id2N = {1};

%% Assembly and Computation
printline(1, 'Assembly and Computation');
while st.next
    s.computeLevel;
end % end stepper iteration

%% Visualization
printline(1, 'Visualization');
s.U.visualize
s.P.visualize
