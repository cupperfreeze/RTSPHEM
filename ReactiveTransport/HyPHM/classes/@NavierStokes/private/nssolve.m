%> @file nssolve.m The raw Navier-Stokes solver.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%>
%>
%> @param  d         Data required by solver [NavierStokes].
%> @param  Uold      Scalar @f$x@f$-velocity of the old time @f$[ \#V+\#E\times 1 ]@f$.
%> @param  Vold      Scalar @f$y@f$-velocity of the old time @f$[ \#V+\#E\times 1 ]@f$.
%> @param  Pold      Scalar pressure of the old time step    @f$[ \#V\times 1 ]@f$.
%> @param  tau       Time step length [scalar].
%> @param  curtime   Evaluation time (of new solution)  [scalar].
%> @param  dataF     Coefficient function for source term in @f$P_2^2@f$ @f$[\#V+\#E\times 2]@f$.
%> @param  dataN     Viscosity, globally constant [scalar].
%> @param  datauD    Dirichlet data for @f$x@f$-velocity @f$[ \#V+\#E\times 1 ]@f$ (def on ALL nodes, but maybe NaN on non-D ones).
%> @param  datavD    Dirichlet data for @f$y@f$-velocity @f$[ \#V+\#E\times 1 ]@f$ (def on ALL nodes, but maybe NaN on non-D ones).
%> @retval  Unew     Scalar @f$x@f$-velocity of the current time @f$[ \#V+\#E\times 1 ]@f$.
%> @retval  Vnew     Scalar @f$y@f$-velocity of the current time @f$[ \#V+\#E\times 1 ]@f$.
%> @retval  Pnew     Scalar pressure of the old time step    @f$[ \#V\times 1 ]@f$.
%> @retval  info       [ struct ]  numIter  [scalar]  number of Newton iterates.
%>
%> @sa HyPHM, NavierStokes, Grid
%> @todo Remove flag isConvection (class NavierStokes not for Stokes any more)

function [U, V, P, info] = nssolve(d, Uold, Vold, tau, dataF, dataN, datauD, datavD, isSlt)

% cf section `MEMORY MANAGEMENT'
persistent storeA; % store stationary (!) matrices such that they can be
persistent storeLocBs; % re-used in following steps
persistent storeLocCs;
persistent storeG;
persistent storeF;
persistent storeNavierStokes; % if problem d does change, delete all stored matrixes

info = [];
wspaces = '                              ';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ASSERTIONS AND ABBREVS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

assert(nargin == 9)
assert(isa(d, 'NavierStokes'))

g = d.grid;
numV = g.numV;
numE = g.numE;
numT = g.numT;
numVE = numV + numE; % abbrevation, often used

assert(length(Uold) == numVE)
assert(length(Vold) == numVE)


%%%%%%%%%%%%%%%%%%%%%%%

%% MEMORY MANAGEMENT %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The persistent variables are still existent when solve() is called from %
% a different instance of NavierStokes.  As remedy, the instance of       %
% NavierStokes is referred to by storeNavierStokes and campared to the    %
% previously used one.  If the previously used does not match all         %
% persistent are deleted.                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(storeNavierStokes)
    storeNavierStokes = d; % very first call of this function
elseif storeNavierStokes ~= d % if d referres to a different problem
    storeNavierStokes = d;
    printline(~isSlt*3, 'Deleting all persistent variables since using a different discretization.')
    storeA = [];
    storeLocBs = [];
    storeLocCs = [];
    storeG = [];
    storeF = [];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  markDirN: [\#V+\#E x 1 logical]

%% DEFINITION OF INDEX SETS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DIRICHLET NODES (for velocity) % for old implementation see F1
markDirE = false(numE, 1); % [\#E x 1] % fetch edge ids acc. to Dirichlet boundary
for k = 1:length(d.id2D)
    markDirE = logical(markDirE+(d.id2D{k} == g.idE));
end

markDirV = false(numV, 1); % [\#V x 1]
for kE = find(markDirE)'
    markDirV(g.V0E(kE, :)) = [true; true]; % intersect. of E_D and E_N is V_D
end
markDirN = [markDirV; markDirE]; % mark all nodes (offset for edge nodes)

% CHECK COMPLETENESS OF BOUNDARY DATA
markNeumE = false(numE, 1); % [\#E x 1] % fetch edge ids acc. to Neumann boundary
for k = 1:length(d.id2N)
    markNeumE = logical(markNeumE+(d.id2N{k} == g.idE));
end
% Boundary has to be well-defined
assert(sum(g.idE ~= 0) == sum(markDirE)+sum(markNeumE), ...
    'HyPHM: Some boundary edges were not assigned to a type.  Please set id2D and/or id2N.');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% UNIQUENESS OF PRESSURE %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~any(markNeumE) && isempty(d.balanceP)
    printline(~isSlt*(-1), ['If there is no Neumann boundary, the pressure balance constraint ', ...
        'should be switched on to ensure uniqueness of pressure.'])
elseif any(markNeumE) && ~isempty(d.balanceP)
    printline(~isSlt*(-1), ['If a pressure balance contraint is defined and there furthermore exists ', ...
        'a Neumann boundary the pressure is overdetermined.'])
end


%%%%%%%%%%%%%%

%% ASSEMBLY %%
%%%%%%%%%%%%%%

tAssembly = tic;
printline(~isSlt*2, 'Assembly of linear terms for a Navier Stokes problem')

%% Assembly of A           <phi_N, phi_N'> (from evaluation term)
if ~isempty(storeA)
    A = storeA; printline(~isSlt*3, 'Using existing A')
else
    A = sparse(numVE, numVE); % quadratic sparse matrix which is build over
    % all indices. Later, the system is solved ona sub index set only.
    if d.isEvolution % only if evolution term is evaluated
        printline(~isSlt*3, 'Assembling stationary A')
        for kT = 1:numT
            nodesOnVerts = g.V0T(kT, :); % vertex numbers on T used for indexing
            nodesOnEdges = g.E0T(kT, :); % edge numbers on T used for indexing
            idxs = [nodesOnVerts, numV + nodesOnEdges]; % using numV as offset
            A(idxs, idxs) = A(idxs, idxs) + locA(d, kT);
        end
    end
    storeA = A; % storing stationary data for next step
end

%% Assembly of B                            B(N',N) = <dx phi_N  dx phi_N'>
printline(~isSlt*3, 'Assembling B')
B = sparse(numVE, numVE); % quadratic sparse matrix which is build over
% all indices. Later, the system is solved on
% a sub index set only.
if ~isempty(storeLocBs)
    locBs = storeLocBs; printline(~isSlt*3, '  Using existing locBs')
else
    printline(~isSlt*3, '  Assembling stationary locBs')
    locBs = cell(numT, 1);
    parfor kT = 1:numT
        locBs{kT} = dataN * locB(d, kT);
    end
    storeLocBs = locBs;
end
printline(~isSlt*3, '  Assembling B via locBs')
for kT = 1:numT
    nodesOnVerts = g.V0T(kT, :); % vertex numbers on T used for indexing
    nodesOnEdges = g.E0T(kT, :); % edge numbers on T used for indexing
    idxs = [nodesOnVerts, numV + nodesOnEdges]; % using numV as offset
    B(idxs, idxs) = B(idxs, idxs) + locBs{kT};
end
clear locBs

%% Assembly of C                            C(N',N) = <dy phi_N  dy phi_N'>
printline(~isSlt*3, 'Assembling C')
C = sparse(numVE, numVE); % quadratic sparse matrix which is build over
% all indices. Later, the system is solved on
% a sub index set only.
if ~isempty(storeLocCs)
    locCs = storeLocCs; printline(~isSlt*3, '  Using existing locCs')
else
    printline(~isSlt*3, '  Assembling stationary locCs')
    locCs = cell(numT, 1);
    parfor kT = 1:numT
        %     warning('may be wrong, check locC')
        locCs{kT} = dataN * locC(d, kT);
    end
    storeLocCs = locCs;
end
printline(~isSlt*3, '  Assembling C via locCs')
for kT = 1:numT
    nodesOnVerts = g.V0T(kT, :); % vertex numbers on T used for indexing
    nodesOnEdges = g.E0T(kT, :); % edge numbers on T used for indexing
    idxs = [nodesOnVerts, numV + nodesOnEdges]; % using numV as offset
    C(idxs, idxs) = C(idxs, idxs) + locCs{kT};
end
clear locCs

%% Assembly of E                            E(N',N) = <dx phi_N  dy phi_N'>
if d.isSymmetric % appears only if D = div (nabla u + (nabla u)')
    printline(~isSlt*3, 'Assembling E')
    E = sparse(numVE, numVE); % quadratic sparse matrix which is build over
    % all indices. Later, the system is solved on
    % a sub index set only.
    for kT = 1:numT
        nodesOnVerts = g.V0T(kT, :); % vertex numbers on T used for indexing
        nodesOnEdges = g.E0T(kT, :); % edge numbers on T used for indexing
        idxs = [nodesOnVerts, numV + nodesOnEdges]; % using numV as offset
        E(idxs, idxs) = E(idxs, idxs) + dataN * locEnum(d, kT);
    end
end

%% Assembly of F (stationary)               F_N'N = <lambda_N', dxphi_N>_kT
if ~isempty(storeF)
    F = storeF; printline(~isSlt*3, 'Using existing F')
else
    printline(~isSlt*3, 'Assembling stationary F')
    F = sparse(numV, numV+numE);
    for kT = 1:numT
        nodesOnVerts = g.V0T(kT, :); % vertex numbers on T used for indexing
        nodesOnEdges = g.E0T(kT, :); % edge numbers on T used for indexing
        idxsU = [nodesOnVerts, numV + nodesOnEdges]; % using numV as offset
        idxsP = nodesOnVerts;
        F(idxsP, idxsU) = F(idxsP, idxsU) + locF(d, kT);
    end
    storeF = F; % storing stationary data for next step
end

%% Assembly of G (stationary)              G_N'N = <lambda_N', dyphi_N'>_kT
if ~isempty(storeG)
    G = storeG; printline(~isSlt*3, 'Using existing G')
else
    printline(~isSlt*3, 'Assembling stationary G')
    G = sparse(numV, numV+numE);
    for kT = 1:numT
        nodesOnVerts = g.V0T(kT, :); % vertex numbers on T used for indexing
        nodesOnEdges = g.E0T(kT, :); % edge numbers on T used for indexing
        idxsU = [nodesOnVerts, numV + nodesOnEdges]; % using numV as offset
        idxsP = nodesOnVerts;

        G(idxsP, idxsU) = ...
            G(idxsP, idxsU) + locG(d, kT);
    end
    storeG = G; % storing stationary data for next step
end

%% Assembly of right hand-side bP
printline(~isSlt*3, 'Assembling bP')
uD = markDirN .* datauD;
vD = markDirN .* datavD; % fetch data on Dirichlet nodes (only)
bP = -F * uD - G * vD;

%% Partial assembly of bU, bV
printline(~isSlt*3, 'Assembling bU and pV (only partially)')
Fxvec = dataF(:, 1);
Fyvec = dataF(:, 2);
stiffF = [1 / 60, -1 / 360, -1 / 360, -1 / 90, 0, 0; ...
    -1 / 360, 1 / 60, -1 / 360, 0, -1 / 90, 0; ...
    -1 / 360, -1 / 360, 1 / 60, 0, 0, -1 / 90; ...
    -1 / 90, 0, 0, 4 / 45, 2 / 45, 2 / 45; ...
    0, -1 / 90, 0, 2 / 45, 4 / 45, 2 / 45; ...
    0, 0, -1 / 90, 2 / 45, 2 / 45, 4 / 45];

bUpart = zeros(numVE, 1); % This is build over all indices.
bVpart = zeros(numVE, 1); % Later, the system is solved on a sub index set only.
for kT = 1:numT
    nodesOnVerts = g.V0T(kT, :); % vertex numbers on T used for indexing
    nodesOnEdges = g.E0T(kT, :); % edge numbers on T used for indexing
    idxs = [nodesOnVerts, numV + nodesOnEdges]; % using numV as offset

    bUpart(idxs) = bUpart(idxs) + det(squeeze(g.A(kT, :, :))) * stiffF * Fxvec(idxs);
    bVpart(idxs) = bVpart(idxs) + det(squeeze(g.A(kT, :, :))) * stiffF * Fyvec(idxs); % old implementation below

    if d.isEvolution
        bUpart(idxs) = bUpart(idxs) + locbUeval(d, kT, tau, Uold);
        bVpart(idxs) = bVpart(idxs) + locbUeval(d, kT, tau, Vold);
    end
end

% write Dirichlet values to solution vector (these entries wont be touched
% by solver)
X = [uD; vD; zeros(numV, 1); 0]; % [U, V, P, dum] where dum is a scalar dummy for the mass balance

printline(~isSlt*3, [wspaces, '...done [%.3f sec]'], toc(tAssembly)) % stop timer assembly linear part


iterErr = inf;
iterCount = 0;
Uiter = Uold;
Viter = Vold;
while iterErr > d.maxRes
    iterCount = iterCount + 1;
    if iterCount > d.maxIter
        error('HyPHM: Maximum number of iteration reached [%d].  Decrease residual [%d].', d.maxIter, d.maxRes)
    end

    %% Assembly of D
    D = sparse(numVE, numVE);
    if d.isConvection % If convection flag is set evaluate additional terms.
        % Otherwise D equals zero.
        for kT = 1:numT
            nodesOnVerts = g.V0T(kT, :); % vertex numbers on T used for indexing
            nodesOnEdges = g.E0T(kT, :); % edge numbers on T used for indexing
            idxs = [nodesOnVerts, numV + nodesOnEdges]; % using numV as offset
            D(idxs, idxs) = D(idxs, idxs) + locD(d, kT, Uiter(idxs), Viter(idxs));
        end
    end

    %% Pressure Balance Contraint               (pressure balance, row vector)
    if ~isempty(d.balanceP)
        mPidx = true;
    else
        mPidx = false;
    end

    % Assembly of mP
    mP = zeros(1, numV);
    if mPidx
        for kT = 1:numT
            idxVs = g.V0T(kT, :);
            mP(idxVs) = mP(idxVs) + g.areaT(kT) / 3;
        end
    end

    % Building right-hand side rhsP
    if mPidx;
        rhsP = sum(g.areaT) * d.balanceP;
    else rhsP = 0;
    end

    %% Build large system of equations
    Z = zeros(1, numVE);
    if d.isSymmetric
        bigA = [1 / tau * A + 2 * B + C + D, E, -F', Z'; ...
            E', 1 / tau * A + B + 2 * C + D, -G', Z'; ...
            F, G, sparse(numV, numV), mP'; ...
            Z, Z, mP, 1]; % additional line for pressure balance
    else
        bigA = [1 / tau * A + B + C + D, sparse(numVE, numVE), -F', Z'; ...
            sparse(numVE, numVE), 1 / tau * A + B + C + D, -G', Z'; ...
            F, G, sparse(numV, numV), mP'; ...
            Z, Z, mP, 1]; % additional line for pressure balance
    end

    %% Clear large variabes if there's no iteration
    if ~d.isConvection
        clear A B C D E G F Fxvec Fyvec
    end
    %   clear store*
    %   whos

    %% Remaining assembly of bU, bV
    bU = bUpart - bigA(1:numVE, 1:2*numVE) * [uD; vD];
    bV = bVpart - bigA(numVE+1:2*numVE, 1:2*numVE) * [uD; vD];

    %% Build large right-hand side
    bigB = [bU; ...
        bV; ...
        bP; ...
        rhsP]; % additional entry for pressure balance

    %% Solving the large system
    % define subindices (solve without Dirichlet indices)
    idx = logical([~markDirN; ~markDirN; ones(numV, 1); mPidx]); % mPidx = false if there's no mass balance

    tSolving = tic; printline(~isSlt*2, 'Solving of the large system for a Navier Stokes problem')
    X(idx) = bigA(idx, idx) \ bigB(idx);
    printline(~isSlt*3, [wspaces, '...done [%.3f sec]'], toc(tSolving)) % stop timer solving

    % extract unknowns from solution vector, set new iterates
    % and compute the error
    U = X(1:numVE);
    V = X(numVE+1:2*numVE);
    P = X(2*numVE+1:2*numVE+numV); % there is the scalar mass balance dummy at the end of X

    % test uniqueness
    dummytol = 1E-13;
    if abs(X(end)) > dummytol
        error('HyPHM: The resulting system is over determined.  Remove constraints!  (...or if dummy ~ tol decrease mesh size) [dummy = %.3e, tol = %.3e]', X(end), dummytol);
    end

    if d.isConvection
        UiterOld = Uiter;
        ViterOld = Viter;
        Uiter = U;
        Viter = V;
        iterErrOld = iterErr;
        iterErr = norm(Uiter-UiterOld) + norm(Viter-ViterOld);
        printline(~isSlt*3, 'Current iteration error: %e (reduction factor: %.2f)', iterErr, iterErrOld/iterErr);
        if iterErr > iterErrOld
            error('HyPHM: Newton''s method for convective term does not converge.  Refine spatial and/or temporal decomposition.')
        end

    else
        iterErr = 0.0; % this ends the while loop
    end


end

printline(~isSlt*3, 'Iteration error: %.3e', iterErr);
printline(~isSlt*3, 'Number of iterations made: %d', iterCount);


%%%%%%%%%%%%%%%%%%%%

%% POSTPROCESSING %%
%%%%%%%%%%%%%%%%%%%%

%% Mean pressure
meanP = 0.0;
for kT = 1:numT
    meanP = meanP + g.areaT(kT) / 3 * sum(P(g.V0T(kT, :)));
end
meanP = meanP / sum(g.areaT);
printline(~isSlt*3, 'Mean pressure of solution is %e', meanP)


return

end


%%%%%%%%%%%%%

%% REMARKS %%
%%%%%%%%%%%%%
%
% The speedup for allocating the sparse matrices is only about 1% => no
% allocation


%%%%%%%%%%%%%%%

%% Footnotes %%
%%%%%%%%%%%%%%%
%
% (F1)
%   %%%  old implementation
% % DIRICHLET NODES (for velocity) % for Neumann nodes see F1
% % mark nodes on vertices
% markDirV = zeros(numV, 1); % [\#V x 1]
% markDirV(unique(g.edgesD)) = 1;
% % mark nodes on edges
% markDirE = zeros(numE, 1); % [\#E x 1]
% if ~isempty(g.edgesD)
%   markDirE(diag(g.V2E(g.edgesD(:, 1), g.edgesD(:, 2)))) ...
%     = ones(size(diag(g.V2E(g.edgesD(:, 1),g.edgesD(:, 2))), 1), 1);
% end
% % mark all nodes (offset for edge nodes)
% markDirN = [markDirV; markDirE];


%   bUpart(idxs) = bUpart(idxs) + locbUsourcenum(d, kT, curtime, 1);
%   bVpart(idxs) = bVpart(idxs) + locbUsourcenum(d, kT, curtime, 2);
