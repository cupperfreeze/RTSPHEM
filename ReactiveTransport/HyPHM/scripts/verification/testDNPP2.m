%> @file testDNPP2.m Convergence test in time for the Darcy-Nernst-Planck-Poisson problem.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
printline(1, 'This is testDNPP2 script.')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters %%
%%%%%%%%%%%%%%%%

% h   = 1/8;  % grid fineness
T = 1; % end time
% tau = 1/8*T;  % step size
tol = 1E-6; % tolerance for fix point iteration
a = 1;
b = 1; % domain lengths
q = 2; % BDFq


levels = 1:4; % 4 is max.
hvec = (1 / 4).^levels;
tauvec = T * (1 / 2).^levels;


isPrintSym = true; % print symbolic definitions?
isVisualizeEx = false; % visualize exact solution?
isVisualizeApprx = false; % visualize approximated solution?
isStore = true; % store the computed errors?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Exact Solution %%
%%%%%%%%%%%%%%%%%%%%
printline(1, 'Building function handles from symbolic solution.')
syms t X Y real % variables
syms alp real % parameter

% Prescribe exact solution
u_sym_1(t, X, Y) = exp(-20*t) * Y;
u_sym_2(t, X, Y) = exp(-20*t) * X;
uN_sym(t, X, Y) = [u_sym_1, u_sym_2] * [2 / a * X - 1; 2 / b * Y - 1];
p_sym(t, X, Y) = exp(-20*t) * cosh(Y) * sin(X);
p_mean_sym(t) = int(int(p_sym(t, X, Y), X, 0, a), Y, 0, b) / (a * b);
c1_sym(t, X, Y) = exp(-10*t) * sinh(X) * cos(Y);
c2_sym(t, X, Y) = exp(-10*t) * cosh(X) * sin(Y);
phi_sym(t, X, Y) = exp(-10*t) * cos(X) * sinh(Y);
% Calculate other components (i.e. no artificial correction term occurs in the resp. equation).
E_sym_1(t, X, Y) = -diff(phi_sym, X);
E_sym_2(t, X, Y) = -diff(phi_sym, Y);

j1_sym_1(t, X, Y) = -diff(c1_sym, X) + (u_sym_1 + E_sym_1) * c1_sym;
j1_sym_2(t, X, Y) = -diff(c1_sym, Y) + (u_sym_2 + E_sym_2) * c1_sym;
j2_sym_1(t, X, Y) = -diff(c2_sym, X) + (u_sym_1 - E_sym_1) * c2_sym;
j2_sym_2(t, X, Y) = -diff(c2_sym, Y) + (u_sym_2 - E_sym_2) * c2_sym;

% Cast as function_handle @(t, x, y, alp).
u_fh_1 = matlabFunction(simplify(u_sym_1));
u_fh_2 = matlabFunction(simplify(u_sym_2));
uN_fh = matlabFunction(simplify(uN_sym));
p_fh = matlabFunction(simplify(p_sym));
p_mean_fh = matlabFunction(simplify(p_mean_sym));
j1_fh_1 = matlabFunction(simplify(j1_sym_1));
j1_fh_2 = matlabFunction(simplify(j1_sym_2));
c1_fh = matlabFunction(simplify(c1_sym));
j2_fh_1 = matlabFunction(simplify(j2_sym_1));
j2_fh_2 = matlabFunction(simplify(j2_sym_2));
c2_fh = matlabFunction(simplify(c2_sym));
E_fh_1 = matlabFunction(simplify(E_sym_1));
E_fh_2 = matlabFunction(simplify(E_sym_2));
phi_fh = matlabFunction(simplify(phi_sym));
Fu_sym_1 = u_sym_1 + diff(p_sym, X) - E_sym_1 * (c1_sym - c2_sym);
Fu_sym_2 = u_sym_2 + diff(p_sym, Y) - E_sym_2 * (c1_sym - c2_sym);
Fu_fh_1 = matlabFunction(simplify(Fu_sym_1));
Fu_fh_2 = matlabFunction(simplify(Fu_sym_2));
Fphi_sym = diff(E_sym_1, X) + diff(E_sym_2, Y) - c1_sym + c2_sym;
Fphi_fh = matlabFunction(simplify(Fphi_sym));
Fc1_sym = diff(c1_sym, t) + diff(j1_sym_1, X) + diff(j1_sym_2, Y);
Fc2_sym = diff(c2_sym, t) + diff(j2_sym_1, X) + diff(j2_sym_2, Y);
Fc1_fh = matlabFunction(simplify(Fc1_sym));
Fc2_fh = matlabFunction(simplify(Fc2_sym));

if isPrintSym
    display(u_fh_1), display(u_fh_2), display(p_fh), ...
        display(j1_fh_1), display(j1_fh_2), display(c1_fh), ...
        display(j2_fh_1), display(j2_fh_2), display(c2_fh), ...
        display(E_fh_1), display(E_fh_2), display(phi_fh);
    display(Fu_fh_1), display(Fu_fh_2), display(Fphi_fh), display(Fc1_fh), display(Fc2_fh);
end


for k = 1:length(levels) % loop over fineness levels


    h = hvec(k);
    g = domainMPP(a, b, h); % [0,1] x [0,1]
    g.print
    tau = tauvec(k);
    st = Stepper(0:tau:T);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Exact and Approximate Solutions on the Grid %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Initialization of exact solution.
    u = Variable(g, st, 'u exact', 'RT0');
    u.setdata(@(t, X, Y) {u_fh_1(t, X, Y); u_fh_2(t, X, Y)});
    p = Variable(g, st, 'p exact', 'P0');
    p.setdata(@(t, X, Y) p_fh(t, X, Y));
    j1 = Variable(g, st, 'j+ exact', 'RT0');
    j1.setdata(@(t, X, Y) {j1_fh_1(t, X, Y); j1_fh_2(t, X, Y)});
    c1 = Variable(g, st, 'c+ exact', 'P0');
    c1.setdata(@(t, X, Y) c1_fh(t, X, Y));
    j2 = Variable(g, st, 'j- exact', 'RT0');
    j2.setdata(@(t, X, Y) {j2_fh_1(t, X, Y); j2_fh_2(t, X, Y)});
    c2 = Variable(g, st, 'c- exact', 'P0');
    c2.setdata(@(t, X, Y) c2_fh(t, X, Y));
    E = Variable(g, st, 'E exact', 'RT0');
    E.setdata(@(t, X, Y) {E_fh_1(t, X, Y); E_fh_2(t, X, Y)});
    phi = Variable(g, st, 'phi exact', 'P0');
    phi.setdata(@(t, X, Y) phi_fh(t, X, Y));
    % Data.
    Fu = Variable(g, st, 'Fu', 'RT0');
    Fu.setdata(@(t, X, Y) {Fu_fh_1(t, X, Y); Fu_fh_2(t, X, Y)}); % correction term in the equation of u
    Fphi = Variable(g, st, 'Fphi', 'P0');
    Fphi.setdata(@(t, X, Y) Fphi_fh(t, X, Y)); % correction term in the equation of phi

    % Approximated solution with start values.
    %   uh    = Variable(g, st, 'uh', 'RT0');    uh.setdata([0,1],   @(t,X,Y)  { u_fh_1(t,X,Y);  u_fh_2(t,X,Y)});
    %   ph    = Variable(g, st, 'ph', 'P0');     ph.setdata([0,1],   @(t,X,Y)  p_fh(t,X,Y));
    %   jh1   = Variable(g, st, 'jh+', 'RT0');   jh1.setdata([0,1],  @(t,X,Y)  {j1_fh_1(t,X,Y); j1_fh_2(t,X,Y)});
    %   ch1   = Variable(g, st, 'ch+', 'P0');    ch1.setdata([0,1],  @(t,X,Y)  c1_fh(t,X,Y));
    %   jh2   = Variable(g, st, 'jh-', 'RT0');   jh2.setdata([0,1],  @(t,X,Y)  {j2_fh_1(t,X,Y); j2_fh_2(t,X,Y)});
    %   ch2   = Variable(g, st, 'ch-', 'P0');    ch2.setdata([0,1],  @(t,X,Y)  c2_fh(t,X,Y));
    %   Eh    = Variable(g, st, 'Eh', 'RT0');    Eh.setdata([0,1],   @(t,X,Y)  { E_fh_1(t,X,Y);  E_fh_2(t,X,Y)});
    %   phih  = Variable(g, st, 'phih', 'P0');   phih.setdata([0,1], @(t,X,Y)  phi_fh(t,X,Y));
    uh = Variable(g, st, 'uh', 'RT0');
    uh.setdata(0, @(t, X, Y) {u_fh_1(t, X, Y); u_fh_2(t, X, Y)});
    ph = Variable(g, st, 'ph', 'P0');
    ph.setdata(0, @(t, X, Y) p_fh(t, X, Y));
    jh1 = Variable(g, st, 'jh+', 'RT0');
    jh1.setdata(0, @(t, X, Y) {j1_fh_1(t, X, Y); j1_fh_2(t, X, Y)});
    ch1 = Variable(g, st, 'ch+', 'P0');
    ch1.setdata(0, @(t, X, Y) c1_fh(t, X, Y));
    jh2 = Variable(g, st, 'jh-', 'RT0');
    jh2.setdata(0, @(t, X, Y) {j2_fh_1(t, X, Y); j2_fh_2(t, X, Y)});
    ch2 = Variable(g, st, 'ch-', 'P0');
    ch2.setdata(0, @(t, X, Y) c2_fh(t, X, Y));
    Eh = Variable(g, st, 'Eh', 'RT0');
    Eh.setdata(0, @(t, X, Y) {E_fh_1(t, X, Y); E_fh_2(t, X, Y)});
    phih = Variable(g, st, 'phih', 'P0');
    phih.setdata(0, @(t, X, Y) phi_fh(t, X, Y));

    %% flow problem
    fp = Transport(g, st, 'Darcy');
    fp.Q = uh;
    fp.U = ph;
    % fp.E is set in the loop.
    %   fp.id2F = {1,2,3,4};
    %   fp.gF.setdata(@(t,X,Y) uN_fh(t,X,Y));
    %  % fp.balanceU is set in loop
    fp.id2D = {1, 2, 3, 4}; % whole boundary is of type Flux
    fp.uD.setdata(@(t, X, Y) p_fh(t, X, Y));

    %% transport problems
    tp1 = Transport(g, st, 'Nernst-Planck 1');
    tp1.Q = jh1;
    tp1.U = ch1;
    % tp1.id2F = {1,3}; % flux
    tp1.id2D = {1, 2, 3, 4}; % Dirichlet
    tp1.A.setdata(@(t, x) 1);
    tp1.F.setdata(@(t, X, Y) Fc1_fh(t, X, Y));
    % tp1.gF.setdata(@(t, x) testDNPP2_g_c(1, t, x));
    tp1.uD.setdata(@(t, X, Y) c1_fh(t, X, Y));

    tp2 = Transport(g, st, 'Nernst-Planck 2');
    tp2.Q = jh2;
    tp2.U = ch2;
    % tp2.id2F = {2,4}; % whole boundary is of type flux
    tp2.id2D = {1, 2, 3, 4}; % whole boundary is of type Dirichlet
    tp2.A.setdata(@(t, x) 1);
    tp2.F.setdata(@(t, X, Y) Fc2_fh(t, X, Y));
    % tp2.gF.setdata(@(t, x) testDNPP2_g_c(2, t, x));
    tp2.uD.setdata(@(t, X, Y) c2_fh(t, X, Y));

    %% electro problem
    ep = Transport(g, st, 'Poisson'); % electro problem
    ep.Q = Eh;
    ep.U = phih;
    % ep.F is set in the loop.
    ep.id2D = {1, 2, 3, 4};
    ep.uD.setdata(@(t, X, Y) phi_fh(t, X, Y));

    st.status
    % %   st.next; %  first time level was initialized above for usage of BDF2
    while st.next % loop over all time steps

        %   ep.balanceU = -2/pi^3*exp(alp*st.curtime)*st.curtime; % <phi> = 2/pi^3 * t;
        ch1data_old = ch1.getdata(st.curstep-1); % initial iterates taken from previous time level
        ch2data_old = ch2.getdata(st.curstep-1);

        itercnt = 0;
        itererr = inf;
        while itererr > tol % fix point interation

            itercnt = itercnt + 1;

            %% electro problem ( set F = c+ - c- )
            ep.F.setdata(st.curstep, ch1data_old-ch2data_old+Fphi.getdata(st.curstep)); % c+ - c-
            ep.computeLevel('silent');

            %% Flow Problem ( set F( eps, c-, c+, qphi) )
            % %       fp.balanceU = p_mean_fh(st.curtime);
            fp.E.setdata(st.curstep, multiplyMatRT0P0(g, eye(2), Eh.getdata(st.curstep), ch1data_old - ch2data_old) ...
                +Fu.getdata(st.curstep));
            fp.computeLevel();

            %% Transport Problems
            tp1.C.setdata(st.curstep, uh.getdata(st.curstep)+Eh.getdata(st.curstep))
            if q == 1 || (q == 2 && st.curstep == 1)
                tp1.computeLevel('silent');
            else % BDF2
                tp1.computeLevel('silent', 'BDF2');
            end
            tp2.C.setdata(st.curstep, uh.getdata(st.curstep)-Eh.getdata(st.curstep))
            if q == 1 || (q == 2 && st.curstep == 1)
                tp2.computeLevel('silent');
            else
                tp2.computeLevel('silent', 'BDF2');
            end

            % IS-error
            ch1data_vold = ch1data_old;
            ch2data_vold = ch2data_old;
            ch1data_old = ch1.getdata(st.curstep);
            ch2data_old = ch2.getdata(st.curstep);

            itererrold = itererr;
            itererr = norm(ch1data_old-ch1data_vold) + norm(ch2data_old-ch2data_vold);

            printline(2, 'Iteration error after iteration %d (t = %.3f, tau = %.3e, h = %.3e)', itercnt, st.curtime, tau, h), printline(3, '%e', itererr)
            if itererr > itererrold
                error('HyPHM: Fixpoint iteration diverges [err: %.1e, errold: %.1e].', itererr, itererrold)
            end
        end % end fix point iteration
        printline(2, 'Timestep and number of iterations')
        printline(3, 'steps: %d,   iterates: %d', st.curstep, itercnt)
    end % end time loop


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Calculation of L2 Error %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    printline(1, 'Computation of L2 errors (tau=1/%d, h=1/%d)', 1/tau, 1/h)

    laststep = st.numsteps; % evaluate L2-error at T

    L2err_u = uh.distance(laststep, @(t, x) [u_fh_1(t, x(1), x(2)); u_fh_2(t, x(1), x(2))], 'L2');
    printline(3, 'L2 error of u: %.8e', L2err_u);
    L2err_p = ph.distance(laststep, @(t, x) p_fh(t, x(1), x(2)), 'L2');
    printline(3, 'L2 error of p: %.8e', L2err_p);
    L2err_j1 = jh1.distance(laststep, @(t, x) [j1_fh_1(t, x(1), x(2)); j1_fh_2(t, x(1), x(2))], 'L2');
    printline(3, 'L2 error of j1: %.8e', L2err_j1);
    L2err_c1 = ch1.distance(laststep, @(t, x) c1_fh(t, x(1), x(2)), 'L2');
    printline(3, 'L2 error of c1: %.8e', L2err_c1);
    L2err_j2 = jh2.distance(laststep, @(t, x) [j2_fh_1(t, x(1), x(2)); j2_fh_2(t, x(1), x(2))], 'L2');
    printline(3, 'L2 error of j2: %.8e', L2err_j2);
    L2err_c2 = ch2.distance(laststep, @(t, x) c2_fh(t, x(1), x(2)), 'L2');
    printline(3, 'L2 error of c2: %.8e', L2err_c2);
    L2err_E = Eh.distance(laststep, @(t, x) [E_fh_1(t, x(1), x(2)); E_fh_2(t, x(1), x(2))], 'L2');
    printline(3, 'L2 error of E: %.8e', L2err_E);
    L2err_phi = phih.distance(laststep, @(t, x) phi_fh(t, x(1), x(2)), 'L2');
    printline(3, 'L2 error of phi: %.8e', L2err_phi);

    printline(3, '%.6e', L2err_u);
    printline(3, '%.6e', L2err_p);
    printline(3, '%.6e', L2err_j1);
    printline(3, '%.6e', L2err_c1);
    printline(3, '%.6e', L2err_j2);
    printline(3, '%.6e', L2err_c2);
    printline(3, '%.6e', L2err_E);
    printline(3, '%.6e', L2err_phi);

    if isStore
        filename = sprintf('testDNPP2_level%d.txt', levels(k));
        dlmwrite(filename, [L2err_u; L2err_p; L2err_j1; L2err_c1; L2err_j2; L2err_c2; L2err_E; L2err_phi], 'precision', '%.8f');
    end

    if isVisualizeEx
        u.visualize, p.visualize, %#ok<UNRCH>
        j1.visualize, c1.visualize,
        j2.visualize, c2.visualize,
        E.visualize, phi.visualize
    end

    if isVisualizeApprx
        uh.visualize, ph.visualize, %#ok<UNRCH>
        jh1.visualize, ch1.visualize,
        jh2.visualize, ch2.visualize,
        Eh.visualize, phih.visualize
    end

end
