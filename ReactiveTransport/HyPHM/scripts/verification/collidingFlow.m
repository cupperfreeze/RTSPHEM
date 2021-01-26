%> @file collidingFlow.m HyPHM script for the test scenario `colliding flow' described in @ref ESW2005.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% disctype = input('Choose NavierStokes [1] or Stokes [2] discretization: ');
disctype = 2;
printline(-1, 'Using the Stokes discretization.  This script can be easily modified to use the N-S discretization.')

%%%%%%%%%%%%%%%%%%%%%%%%%

%% Analytical Solution %%
%%%%%%%%%%%%%%%%%%%%%%%%%

uAnal = @(t, x) 20 * x(1) * x(2)^3;
vAnal = @(t, x) 5 * x(1)^4 - 5 * x(2)^4;
pAnal = @(t, x) 60 * x(1)^2 * x(2) - 20 * x(2)^3;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization of Timer %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

st = Stepper([0, 1]); % one single step from 0 to 1
st.print
% st.status


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization of Grids %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hmax = input('Define h: ');
g = domainRectangle(-1, 1, -1, 1, hmax);
% g.visualize
g.print


%%%%%%%%%%%%%%%%%%%%%%

%% Solution Vectors %%
%%%%%%%%%%%%%%%%%%%%%%

pres = Variable(g, st, 'pressure', 'P1');
velo = Variable(g, st, 'velocity', 'P2P2');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization of Problem Data %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if disctype == 1
    cf = NavierStokes(g, st); printline(3, 'Choosing NavierStokes discretization.')

    % Flags
    cf.isConvection = false; % convective term?
    cf.isEvolution = false; % evolution term?
    cf.isSymmetric = false; % symmetric diff tensor?
elseif disctype == 2
    cf = Stokes(g, st, 'Colliding Flow');
    printline(3, 'Choosing Stokes discretization.')
    cf.isEvolution = false; % evolution term?
end

% Boundary
cf.id2D = {1, 2, 3, 4};

% Mean Pressure Contraint
cf.balanceP = 0; % mean pressure equal zero?

% Parameters and Coefficient Functions
visc = Variable(g, st, 'viscosity', 'C');
visc.setdata(1);
cf.N = visc;
source = Variable(g, st, 'source', 'P2P2');
source.setdata(@(t, x) zeros(2, 1));
cf.F = source;


bdryflow = Variable(g, st, 'bdryflow', 'P2P2');
bdryflow.setdata(@(t, x) [uAnal(t, x); vAnal(t, x)]);
cf.uD = bdryflow;
cf.print


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization of Initial Solution %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initializaton not needed here, but still quoted for reference
pres.setdata(0, zeros(g.numV, 1));
velo.setdata(0, zeros(g.numV + g.numE, 2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Assembly and Computation %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

printline(1, 'Assembly and Computation');
timer_computation = clock;

while st.next

    cf.computeLevel;

end % end stepper iteration

printline(2, 'Total assembly and computation time');
printline(3, '%.2f sec', etime(clock, timer_computation));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculation of L2 Error %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% split velocity in x and y-components, since P2 has a distance-method
u = Variable(g, st, '', 'P2');
v = Variable(g, st, '', 'P2');
velodata = cf.U.getdata(st.numsteps);
u.setdata(1, velodata(:, 1));
v.setdata(1, velodata(:, 2));

L2err_u = u.distance(1, uAnal, 'L2');
L2err_v = v.distance(1, vAnal, 'L2');
L2err_p = cf.P.distance(1, pAnal, 'L2');

printline(3, 'L2 error of u: %.3e', L2err_u);
printline(3, 'L2 error of v: %.3e', L2err_v);
printline(3, 'L2 error of p: %.3e', L2err_p);


%%%%%%%%%%%%%%%%%%%

%% Visualization %%
%%%%%%%%%%%%%%%%%%%

printline(1, 'Visualization');

cf.U.visualize
cf.P.visualize
