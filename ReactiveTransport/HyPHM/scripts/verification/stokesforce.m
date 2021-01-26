%> @file stokesforce.m HyPHM script for to test the right-hand side in Stokes equation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%

%% Analytical Solution %%
%%%%%%%%%%%%%%%%%%%%%%%%%

% collidingFlow is stationary, time is only set for reference purposses
uAnal = @(x) -cos(pi*x(1)) * sin(pi*x(2));
vAnal = @(x) sin(pi*x(1)) * cos(pi*x(2));
pAnal = @(x) 0;
% the inner force is then given by
fAnal = @(t, x) 2 * pi^2 * [uAnal(x); vAnal(x)];
% the Dirichlet data is obtained by sampling the solution on the boundary
uDAnal = @(t, x) [uAnal(x); vAnal(x)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization of Timer %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

st = Stepper([0, 1]); % one single step from 0 to 1
st.print


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization of Grids %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hmax = input('Define h: ');
g = domainRectangle(0, 1, 0, 1, hmax);
% g.visualize
g.print


%%%%%%%%%%%%%%%%%%%%%%

%% Solution Vectors %%
%%%%%%%%%%%%%%%%%%%%%%

force = Variable(g, st, 'force', 'P2P2');
force.setdata(fAnal)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization of Problem Data %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sf = Stokes(g, st, 'stokesforce');
sf.isEvolution = false; % evolution term?

% Boundary
sf.id2D = {1, 2, 3, 4};

% Mean Pressure Contraint
sf.balanceP = 0; % mean pressure equal zero?

% Parameters and Coefficient Functions
force = Variable(g, st, 'force', 'P2P2');
force.setdata(fAnal);
sf.F = force;

bdryforces = Variable(g, st, '', 'P2P2');
bdryforces.setdata(uDAnal);
sf.uD = bdryforces;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Assembly and Computation %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

printline(1, 'Assembly and Computation');
timer_computation = clock;

while st.next, sf.computeLevel, end % end stepper iteration

printline(2, 'Total assembly and computation time');
printline(3, '%.2f sec', etime(clock, timer_computation));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculation of L2 Error %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% split velocity in x and y-components, since P2 has a distance-method
u = Variable(g, st, '', 'P2');
v = Variable(g, st, '', 'P2');
velodata = sf.U.getdata(st.numsteps);
u.setdata(1, velodata(:, 1));
v.setdata(1, velodata(:, 2));

L2err_u = u.distance(1, @(t, x) uAnal(x), 'L2');
L2err_v = v.distance(1, @(t, x) vAnal(x), 'L2');
L2err_p = sf.P.distance(1, @(t, x) pAnal(x), 'L2');

printline(3, 'L2 error of u: %.3e', L2err_u);
printline(3, 'L2 error of v: %.3e', L2err_v);
printline(3, 'L2 error of p: %.3e', L2err_p);

printline(3, '\nDOF = %d, h = %f.', g.numE*2+g.numV*3, hmax)
