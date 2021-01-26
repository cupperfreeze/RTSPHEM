%> @file stdMiscibleDisplacement.m Standard script for solving the equation for miscible displacement.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters
muS = 0.001; % viscosity of solvent
muR = 1; % viscosity of resident fluid
mu = @(c) ((c * muS).^(-1 / 4) + ((1 - c) * muR).^(-1 / 4)).^(-4); % viscosity

Dm = 1E-5; % molecular diffusivity
alphaL = 1; % longitudinal dispersion coefficient
alphaT = 1; % transversal dispersion coefficient

hmax = 5E-2; % mesh fineness
tau = 1E-3; % time step size
Tend = 1000 * tau;

%%

%% Grid
printline(1, 'Grid Generation')
g = domainRectangle(0, 1, 0, 1, hmax);
g = FoldedGrid(g, 1);
g.print
% g.visualize('idE')

%% More Parameters
phi = rand(g.numT, 1) * 0.5 + 0.1; % porosity
K = 1E-1 * phi; % permeability of the medium

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

printline(1, 'Pre-processing')

%% Stepper
st = Stepper(0:tau:Tend);
st.print

%% Unknowns
uh = Variable(g, st, 'water flux', 'RT0');
uh.setdata(0, @(t, x) [0.0; 0.0]);
ph = Variable(g, st, 'pressure', 'P0');
ph.setdata(0, @(t, x) 0.0);
jh = Variable(g, st, 'mass flux', 'RT0');
jh.setdata(0, @(t, x) [0.0; 0.0]);
ch = Variable(g, st, 'concentration', 'P0');
ch.setdata(0, @(t, x) 0.02+0.90*((x(1) <= 0.2) || (norm(x - [0.2; 0.5]) <= 0.2)));

%% Darcy Problem
darcy = Transport(g, st, 'Darcy problem');
% unknowns
darcy.Q = uh;
darcy.U = ph;
% coefficients
darcy.D = Variable(g, st, 'permeability', 'P0P0P0P0');
darcy.D.setdata(0, zeros(g.numT, 4)); % computeLevel requires the diffusion at the old time step
% boundary data
darcy.id2F = {4};
darcy.id2D = {2};
darcy.gF.setdata(@(t, x) -0.01);
darcy.uD.setdata(@(t, x) 0.0);

%% Transport Problem
trans = Transport(g, st, 'Transport problem');
% unknowns
trans.Q = jh;
trans.U = ch;
% coefficients
trans.A.setdata(phi);
trans.D = Variable(g, st, 'diffusion', 'P0P0P0P0');
trans.D.setdata(0, zeros(g.numT, 4)); % computeLevel requires the diffusion at the old time step
% boundary data
trans.id2D = {4};
trans.id2N = {2};
trans.uD.setdata(@(t, x) (x(1) == 0)*0.95);
% flags
trans.isUpwind = 'alt';

%% Assembly and Computation
printline(1, 'Assembly and Computations')
st.status
while st.next


    %   (0.5*sin(10*x(1))*sin(10*x(2)) + 0.51)
    Kdata = [1 ./ mu(ch.getdata(st.curstep - 1)) .* K, zeros(g.numT, 2), 1 ./ mu(ch.getdata(st.curstep - 1)) .* K];
    max(Kdata(:))
    min(Kdata(:))
    darcy.D.setdata(st.curstep, Kdata);
    darcy.computeLevel;

    udata = uh.getdata(st.curstep, 'P0P0'); % [#T, 2]
    Ddata = zeros(g.numT, 4);
    for kT = 1:g.numT
        uh_T = udata(kT, :)';
        D_T = (alphaT * norm(uh_T) + Dm) * eye(2) + (alphaL - alphaT) * (uh_T * uh_T') / norm(uh_T);
        Ddata(kT, :) = D_T(:);
    end
    trans.D.setdata(st.curstep, Ddata);
    trans.C.setdata(st.curstep, uh.getdata(st.curstep));
    trans.computeLevel;
end

%% Visualization
printline(1, 'Visualization')
uh.visualize % water flux
ph.visualize % pressure head
ch.visualize % concentration