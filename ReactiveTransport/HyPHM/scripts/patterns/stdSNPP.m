%> @file stdSNPP.m Stokes&ndash;Nernst&ndash;Planck&ndash;Poisson system with one species.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%>
%> Solves the fully-coupled Stokes–Nernst–Planck–Poisson system on the pore scale
%> including physical quantities/dimensions.
%>


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Triangulation %%
%%%%%%%%%%%%%%%%%%%
len = 1E-5; % domain size
h = .03 * len; % mesh fineness
tau = .005; % step size
T = .1; % end time
tol = 1E-6; % tolerance for fixed-point iteration


rotang = 0; % rotation angle of ellipsis
gd = [3, 4; ... % geometry description
    4, len / 2; ...
    0, len / 2; ...
    len, len / 6; ...
    len, len / 6; ...
    0, rotang; ...
    0, 0; ...
    0, 0; ...
    len, 0; ...
    len, 0];
sf = 'square-ellipse'; % set formula
ns = names2ns('square', 'ellipse'); % name space
[p, e, t] = initmesh(decsg(gd, sf, ns), 'Hmax', h);
g = Grid(p, e, t);
% g.visualize('idE')
g = FoldedGrid(g, 'semifold');
st = Stepper(0:tau:T);
% boundary identification
dOin = {3}; % inflow boundary
dOout = {1}; % outflow boundary
dOsurf = {2, 4}; % surface/impermeable in \partial\Omega
gam = num2cell(5:1E3); % interior boundary Gamma


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Constants and Pseudoconstants at 20°C %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C_NA = 6.022E+23; % Avogadro number [mol^{-1}]
C_kB = 1.381E-23; % Boltzmann constant [J K^{-1}]
C_e = 1.602E-19; % elementary charge [C]
C_F = 9.648E+4; % Faraday constant [C mol^{-1}]
C_R = 8.314E+0; % gas constant [J K^{-1} mol^{-1}]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C_eps = 5.553E-8; % permittivity (of water) [C V^{-1} m^{-1}]
C_mu = 1.002E-3; % dynamic viscosity (of water) [Pa s]
C_nu = 1.004E-6; % kinematic viscosity (of water) [m^2 s^{-1}]
C_rho = 9.982E+2; % water density [kg\,m^{-3}]
C_T = 293; % temperature [K]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Data %%
%%%%%%%%%%
u_in = @(t, x) 1E-3 * (x(1) == 0) * [1; 0]; % [m s^{-1}]  Y. Liu ``Transport of Cryptosporidium parvum...'' used 2.9E-3 m/s
c_in = @(t, x) 1E-4; % [mol m^{-3}]
sigma = @(t, x) 0; % surface charge density [C m^{-2}]
phi_D = @(t, x) 0; % surface potential [V]
v = 4.105086941636337E-7; % electrical mobility [mol s kg^{-1}]
diff  = C_R*C_T*v; % diffusivity [m^2 s^{-1}]
z = 1; % valence of the species [-]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Unknowns %%
%%%%%%%%%%%%%%
u = Variable(g, st, 'u', 'P2P2');
u.setdata(0, @(t, x) [0; 0]);
p = Variable(g, st, 'p', 'P1');
p.setdata(0, @(t, x) 0);
j   = Variable(g, st, 'j',   'RT0');    j.setdata(0,  @(t, x) [0;0]);
c = Variable(g, st, 'c', 'P0');
c.setdata(0, @(t, x) 0);
E = Variable(g, st, 'E', 'RT0');
E.setdata(0, @(t, x) [0; 0]);
phi = Variable(g, st, 'phi', 'P0');
phi.setdata(0, @(t, x) 0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Stokes Problem %%
%%%%%%%%%%%%%%%%%%%%
p_S = Stokes(g, st, 'Stokes problem');
p_S.U = u;
p_S.P = p;
p_S.id2D = [dOin, dOsurf, gam]; % inflow left, no-flow at obstacles
p_S.uD = Variable(g, st, '', 'P2P2');
p_S.uD.setdata(u_in);
p_S.id2N = dOout; % outflow right
p_S.N = Variable(g, st, '', 'C');
p_S.N.setdata(C_eps);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Nernst–Planck Problem %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
p_NP = Transport(g, st, 'Nernst–Planck problem');
p_NP.Q = j;
p_NP.U = c;
p_NP.id2F = [dOsurf, gam]; % top & bottom are impermeable and so are the interiors
p_NP.gF = Variable(g, st, '', 'P0E');
p_NP.gF.setdata(@(t, x) 0);
p_NP.id2D = dOin;
p_NP.uD = Variable(g, st, '', 'P0E');
p_NP.uD.setdata(c_in);
p_NP.id2N = dOout; % right bdry is outflow
p_NP.A.setdata(@(t, x) 1);
p_NP.D.setdata(eye(2)*C_R*C_T*v);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Poisson Problem %%
%%%%%%%%%%%%%%%%%%%%%
p_P = Transport(g, st, 'Poisson problem'); % electro problem
p_P.Q = E;
p_P.U = phi;
p_P.id2F = [dOin, dOout, gam];
p_P.gF.setdata(@(t, x) sigma(t, x)/C_eps);
p_P.id2D = dOsurf;
p_P.uD.setdata(phi_D);
p_P.D.setdata(eye(2));


globitercnt = 0;
st.status % show status bar
while st.next % loop over all time stp_Ps

    cdata_iter = c.getdata(st.curstep-1); % initial iterates taken from previous time level

    itercnt = 0;
    itererr = inf;
    while itererr > tol % fix point interation

        itercnt = itercnt + 1;

        %% Poisson Problem
        p_P.F.setdata(st.curstep, (C_F * z / C_eps)*cdata_iter);
        p_P.computeLevel('silent')

        %% Flow Problem ( set F( p_Ps, c-, c+, qphi) )
        %     Ehdata_P0P0 = RT0.RT0toP0P0slice(g, E.getdata(st.curstep)); % [#T x 2]
        %     p_S.F.setdata(st.curstep,  (C_F*z) * P0P0.P0P0toP2P2slice(g, [cdata_iter.*Ehdata_P0P0(:,1), cdata_iter.*Ehdata_P0P0(:,2)]) ) % F*z/eps*E*c
        p_S.computeLevel('silent')

        %% Tranp_Sort Problem
        %     p_NP.C.setdata(st.curstep, u.getdata(st.curstep, 'RT0') + ((diff*z*C_e)/(C_kB*C_T))*E.getdata(st.curstep)) % u + D*z*e/kB/T*E
        p_NP.C.setdata(st.curstep, u.getdata(st.curstep, 'RT0')) % u + D*z*e/kB/T*E
        p_NP.computeLevel()

        % IS-error
        cdata_iter_old = cdata_iter;
        cdata_iter = c.getdata(st.curstep);

        itererrold = itererr;
        itererr = norm(cdata_iter-cdata_iter_old);

        printline(2, 'Current iteration error'), printline(3, '%e', itererr)
        %  itererr = 0
        if itererr > itererrold
            error('HyPHM: Fixpoint iteration diverges at step %d [err: %.8e, errold: %.8e].', itercnt, itererr, itererrold)
        end
    end % end fix point iteration
    printline(2, 'Timestp_P and number of iterations')
    printline(3, 'steps: %d,   iterates: %d', st.curstep, itercnt)
    globitercnt = globitercnt + itercnt;
end % end time loop

printline(3, 'Total number of iterations: %d', globitercnt)
printline(3, 'Average number of iterations per step: %.3f', globitercnt/st.numsteps)

u.visualize
% p.visualize
% j.visualize
c.visualize
E.visualize
% phi.visualize
