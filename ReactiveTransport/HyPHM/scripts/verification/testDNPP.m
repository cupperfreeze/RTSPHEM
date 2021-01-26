%> @file testDNPP.m Script for the Darcy-Nernst-Planck-Poisson problem.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = 1 / 16; % grid fineness
tau = 1 / 8; % step size
T = 1; % end time
tol = 1E-14; % tolerance for fix point iteration

isVisualizeEx = false; % visualize exact solution?
isVisualizeApprx = false; % visualize approximated solution?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% grid
g = domainMPP(1, 1, h); % [0,1] x [0,1]
% g.visualize('idE')
g.print

%% timer
st = Stepper(0:tau:T);

%% exact solutions
u = Variable(g, st, 'u exact', 'RT0');
u.setdata(@testDNPP_u);
p = Variable(g, st, 'p exact', 'P0');
p.setdata(@testDNPP_p);
force = Variable(g, st, 'force', 'RT0');
force.setdata(@testDNPP_force);
q1 = Variable(g, st, 'q+ exact', 'RT0');
q1.setdata(@(t, x) testDNPP_q(1, t, x));
c1 = Variable(g, st, 'c+ exact', 'P0');
c1.setdata(@(t, x) testDNPP_c(1, t, x));
q2 = Variable(g, st, 'q- exact', 'RT0');
q2.setdata(@(t, x) testDNPP_q(2, t, x));
c2 = Variable(g, st, 'c- exact', 'P0');
c2.setdata(@(t, x) testDNPP_c(2, t, x));
E = Variable(g, st, 'E exact', 'RT0');
E.setdata(@testDNPP_E);
phi = Variable(g, st, 'phi exact', 'P0');
phi.setdata(@testDNPP_phi);

%% approximated solution + initialization
uh = Variable(g, st, 'uh', 'RT0');
uh.setdata(0, @testDNPP_u);
ph = Variable(g, st, 'ph', 'P0');
ph.setdata(0, @testDNPP_p);
qh1 = Variable(g, st, 'qh+', 'RT0');
qh1.setdata(0, @(t, x) testDNPP_q(1, t, x));
ch1 = Variable(g, st, 'ch+', 'P0');
ch1.setdata(0, @(t, x) testDNPP_c(1, t, x));
qh2 = Variable(g, st, 'qh-', 'RT0');
qh2.setdata(0, @(t, x) testDNPP_q(2, t, x));
ch2 = Variable(g, st, 'ch-', 'P0');
ch2.setdata(0, @(t, x) testDNPP_c(2, t, x));
Eh = Variable(g, st, 'Eh', 'RT0');
Eh.setdata(0, @testDNPP_E);
phih = Variable(g, st, 'phih', 'P0');
phih.setdata(0, @testDNPP_phi);

%% electro problem
ep = Transport(g, st, 'Poisson'); % electro problem
ep.Q = Eh;
ep.U = phih;
ep.id2F = {1, 2, 3, 4}; % whole boundary is of type flux/Neumann
ep.gF.setdata(@testDNPP_g_phi);

%% flow problem
fp = Transport(g, st, 'Darcy');
fp.Q = uh;
fp.U = ph;
fp.balanceU = 0;
fp.id2F = {1, 2, 3, 4}; % whole boundary is of type Flux
fp.gF.setdata(@testDNPP_uN);

%% transport problems
tp1 = Transport(g, st, 'Nernst-Planck 1');
tp1.Q = qh1;
tp1.U = ch1;
tp1.id2F = {1, 3}; % flux
tp1.id2D = {2, 4}; % Dirichlet
tp1.A.setdata(@(t, x) 1);
tp1.F.setdata(@(t, x) testDNPP_s(1, t, x));
tp1.gF.setdata(@(t, x) testDNPP_g_c(1, t, x));
tp1.uD.setdata(@(t, x) testDNPP_c(1, t, x));

tp2 = Transport(g, st, 'Nernst-Planck 2');
tp2.Q = qh2;
tp2.U = ch2;
tp2.id2F = {2, 4}; % whole boundary is of type flux
tp2.id2D = {1, 3}; % whole boundary is of type Dirichlet
tp2.A.setdata(@(t, x) 1);
tp2.F.setdata(@(t, x) testDNPP_s(2, t, x));
tp2.gF.setdata(@(t, x) testDNPP_g_c(2, t, x));
tp2.uD.setdata(@(t, x) testDNPP_c(2, t, x));

st.status
while st.next % loop over all time steps

    ep.balanceU = -2 / pi^3 * st.curtime; % <phi> = 2/pi^3 * t;

    ch1data_old = ch1.getdata(st.curstep-1); % initial iterates taken from previous time level
    ch2data_old = ch2.getdata(st.curstep-1);

    itercnt = 0;
    itererr = inf;
    while itererr > tol % fix point interation

        itercnt = itercnt + 1;

        %% electro problem ( set F = c+ - c- )
        ep.F.setdata(st.curstep, ch1data_old-ch2data_old); % c+ - c-
        ep.computeLevel('silent');

        %% Flow Problem ( set F( eps, c-, c+, qphi) )
        fp.E.setdata(st.curstep, multiplyMatRT0P0(g, eye(2), Eh.getdata(st.curstep), ch1data_old - ch2data_old) ...
            +force.getdata(st.curstep));
        fp.computeLevel('silent');

        %% Transport Problems
        tp1.C.setdata(st.curstep, uh.getdata(st.curstep)+Eh.getdata(st.curstep))
        tp1.computeLevel('silent');

        tp2.C.setdata(st.curstep, uh.getdata(st.curstep)-Eh.getdata(st.curstep))
        tp2.computeLevel('silent');

        % IS-error
        ch1data_vold = ch1data_old;
        ch2data_vold = ch2data_old;
        ch1data_old = ch1.getdata(st.curstep);
        ch2data_old = ch2.getdata(st.curstep);

        itererrold = itererr;
        itererr = norm(ch1data_old-ch1data_vold) + norm(ch2data_old-ch2data_vold);

        printline(2, 'Iteration error after iteration %d (tau = %.3e, h = %.3e)', itercnt, tau, h), printline(3, '%e', itererr)
        if itererr > itererrold
            error('HyPHM: Fixpoint iteration diverges [err: %.1e, errold: %.1e].', itererr, itererrold)
        end
    end % end fix point iteration
    printline(2, 'Timestep and number of iterations')
    printline(3, 'steps: %d,   iterates: %d', st.curstep, itercnt)
end % end time loop

if isVisualizeEx
    u.visualize, p.visualize, q1.visualize, c1.visualize, q2.visualize, c2.visualize, E.visualize, phi.visualize
end

if isVisualizeApprx
    uh.visualize, ph.visualize, qh1.visualize, ch1.visualize, qh2.visualize, ch2.visualize, Eh.visualize, phih.visualize
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculation of L2 Error %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

printline(1, 'Computation of L2 errors')

laststep = st.numsteps; % evaluate L2-error at T

L2err_u = uh.distance(laststep, @testDNPP_u, 'L2');
printline(3, 'L2 error of u: %.8e', L2err_u);
L2err_p = ph.distance(laststep, @testDNPP_p, 'L2');
printline(3, 'L2 error of p: %.8e', L2err_p);
L2err_q1 = qh1.distance(laststep, @(t, x) testDNPP_q(1, t, x), 'L2');
printline(3, 'L2 error of q1: %.8e', L2err_q1);
L2err_c1 = ch1.distance(laststep, @(t, x) testDNPP_c(1, t, x), 'L2');
printline(3, 'L2 error of c1: %.8e', L2err_c1);
L2err_q2 = qh2.distance(laststep, @(t, x) testDNPP_q(2, t, x), 'L2');
printline(3, 'L2 error of q2: %.8e', L2err_q2);
L2err_c2 = ch2.distance(laststep, @(t, x) testDNPP_c(2, t, x), 'L2');
printline(3, 'L2 error of c2: %.8e', L2err_c2);
L2err_E = Eh.distance(laststep, @testDNPP_E, 'L2');
printline(3, 'L2 error of E: %.8e', L2err_E);
L2err_phi = phih.distance(laststep, @testDNPP_phi, 'L2');
printline(3, 'L2 error of phi: %.8e', L2err_phi);

printline(3, '%.6e', L2err_u);
printline(3, '%.6e', L2err_p);
printline(3, '%.6e', L2err_q1);
printline(3, '%.6e', L2err_c1);
printline(3, '%.6e', L2err_q2);
printline(3, '%.6e', L2err_c2);
printline(3, '%.6e', L2err_E);
printline(3, '%.6e', L2err_phi);


g.print
display(h)