%> @file RBPA2008.m HyPHM script for the test scenario appearing in @ref RBPA2008.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> Related files are in scripts/RBPA2008
%>
%> Dependencies:
%>   - RBPA2008_u.m
%>   - RBPA2008_q.m
%>   - RBPA2008_F.m
%>   - errorRBPA2008.m


addpath('./scripts/scripts-discerrors/RBPA2008') % add folder to path where scenario stuff is located

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization of Grids %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hmax = input('Define h: ');
g = Grid('RBPA2008.geo', hmax);
% g.visualize
g.print


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization of Timer %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

st = Stepper(0:.1:1);
st.print
st.status


%%%%%%%%%%%%%%%%%%%%%%

%% Solution Vectors %%
%%%%%%%%%%%%%%%%%%%%%%

conc1 = Unknown(g, st, 'concentration acceptor', 'P0');
mflux1 = Unknown(g, st, 'mass flux acceptor', 'RT0');
conc2 = Unknown(g, st, 'concentration donator', 'P0');
mflux2 = Unknown(g, st, 'mass flux donator', 'RT0');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization of Problem Data %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kForw = 1.0;
kBackw = 0.0;
stoichMat = [2; 1];
d = ConvDiff(g, st, stoichMat, kForw, kBackw);
d.id2D = {1, 2, 3, 4}; % Boundary

d.A = @(t, x, specNo) 1.0;
d.B = @(t, x, specNo) 0.0;
d.C = @(t, x, specNo) [0.0; -1.0]; % constant vertical flow (no coupling)
d.D = @(t, x, specNo) 1E-1 * eye(2);
d.E = @(t, x, specNo) [0.0; 0.0];
d.F = @RBPA2008_F;
d.uD = @RBPA2008_u;

d.maxRes = 1e-15;
d.maxIter = 1E3;

d.print


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization of Initial Solution %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

transInit = @RBPA2008_u;

conc1init = zeros(g.numT, 1);
conc2init = zeros(g.numT, 1);
for kT = 1:g.numT
    conc1init(kT) = transInit(0.0, g.baryT(kT, :)', 1);
    conc2init(kT) = transInit(0.0, g.baryT(kT, :)', 2);
end

conc1.setdata(0, conc1init);
conc2.setdata(0, conc2init);
mflux1.setdata(0);
mflux2.setdata(0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Assembly and Computation %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

printline(1, 'Assembly and Computation');
tic;

while st.next;
    d.computeLevel({mflux1, mflux2}, {conc1, conc2});
end % end stepper iteration

printline(2, 'Total assembly and computation time')
printline(3, '%.2f sec', toc)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculation of L2 Error %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

errorRBPA2008


%%%%%%%%%%%%%%%%%%%

%% Visualization %%
%%%%%%%%%%%%%%%%%%%

printline(1, 'Visualization')

conc1.visualize
conc2.visualize
mflux1.visualize
mflux2.visualize
