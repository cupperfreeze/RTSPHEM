%> @file BC2005.m HyPHM script for the test scenario appearing in @ref BC2005.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%>
%> Related files are in
%>   sripts/BC2005
%>
%> Dependencies:
%>   - BC2005_gF.m
%>   - BC2005_q.m
%>   - BC2005_q_disp.m  % display solution via quiver
%>   - BC2005_u.m
%>   - errorBC2005.m
%>

addpath('./scripts/scripts-discerrors/BC2005') % add folder to path where scenario stuff is located

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization of Timer %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

st = Stepper([0, 1]); % one single step from zero to one
st.print


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization of Grids %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calling the mesh generator
printline(2, 'Please define max. meshsize <= 1')
hmax = input('');

% Create new Grid instance
g = Grid('BC2005.geo', hmax);

% Print information
g.print
% g.visualize('numE', 'numT')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization of Problem Data %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = Transport(g, st, 'transport');
t.id2D = {2, 3};
t.id2F = {1, 4, 5, 6};

t.A = Variable(g, st, 'A', 'P0');
t.A.setdata(@(t, x) 0.0);
t.B = Variable(g, st, 'B', 'P0');
t.B.setdata(@(t, x) 0.0);
t.C = Variable(g, st, 'C', 'RT0');
t.C.setdata(zeros(g.numE, 1));
t.D = Variable(g, st, 'D', 'CCCC');
t.D.setdata(eye(2));
t.E = Variable(g, st, 'E', 'RT0');
t.E.setdata(zeros(g.numE, 1));
t.F = Variable(g, st, 'F', 'P0');
t.F.setdata(@(t, x) 0.0);
t.uD = Variable(g, st, 'uD', 'P0E');
t.uD.setdata(@(t, x) 0.0);

gFdata = zeros(g.numE, 1); % Neumann is identical to t.Q here.  Thus t.Q is
% used since only homog. Neum. was implemented.
for kE = 1:g.numE
    gFdata(kE) = BC2005_gF(g.baryE(kE, :)', g.nuE(kE, :)');
end
t.gF = Variable(g, st, 'gF', 'P0E');
t.gF.setdata(gFdata);

t.print

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Assembly and Computation %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

printline(1, 'Assembly and Computation')
timer_computation = clock;

while st.next

    t.computeLevel();

end % end stepper iteration

printline(2, 'Total assembly and computation time')
printline(3, '%.2f sec', etime(clock, timer_computation))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculation of L2 Error %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

err_u_L2 = t.U.distance(1, @(t, x) BC2005_u(1, x, 1), 'L2');
err_u_P0 = t.U.distance(1, @(t, x) BC2005_u(1, x, 1), 'P0');
err_q_L2 = t.Q.distance(1, @(t, x) BC2005_q(1, x, 1), 'L2');
err_q_RT0 = t.Q.distance(1, @(t, x) BC2005_q(1, x, 1), 'RT0');

printline(3, 'L2 error of %s at t = %.1f: %.3e', t.U.name, 1, err_u_L2);
printline(3, 'P0 error of %s at t = %.1f: %.3e', t.U.name, 1, err_u_P0);
printline(3, 'L2 error of %s at t = %.1f: %.3e', t.Q.name, 1, err_q_L2);
printline(3, 'RT0 error of %s at t = %.1f: %.3e', t.Q.name, 1, err_q_RT0);


%%%%%%%%%%%%%%%%%%%

%% Visualization %%
%%%%%%%%%%%%%%%%%%%

printline(1, 'Visualization')
t.U.visualize
t.Q.visualize
