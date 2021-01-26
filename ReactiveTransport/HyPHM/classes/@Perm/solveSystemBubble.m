%> @file Perm/solveSystem.m Solve the linear system for the Stokes problem if Discretization = Bubble.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @brief Solve the linear system for the Stokes problem.
%>
%> @param  this      [Perm].
%> @param  Uold      Scalar @f$x@f$-velocity of the old time @f$[ \#V+\#E\times 1 ]@f$.
%> @param  Vold      Scalar @f$y@f$-velocity of the old time @f$[ \#V+\#E\times 1 ]@f$.
%> @param  tau       Time step length [scalar].
%> @param  dataN     Viscosity, globally constant [scalar].
%> @param  datauD    Dirichlet data for @f$x@f$-velocity @f$[ \#V+\#E\times 1 ]@f$ (def on ALL nodes, but maybe NaN on non-D ones).
%> @param  datavD    Dirichlet data for @f$y@f$-velocity @f$[ \#V+\#E\times 1 ]@f$ (def on ALL nodes, but maybe NaN on non-D ones).
%> @retval  U        Scalar @f$x@f$-velocity of the current time @f$[ \#V+\#T\times 1 ]@f$ for force term e_1.
%> @retval  V        Scalar @f$y@f$-velocity of the current time @f$[ \#V+\#T\times 1 ]@f$ for force term e_1.
%> @retval  U2       Scalar @f$x@f$-velocity of the current time @f$[ \#V+\#T\times 1 ]@f$ for force term e_2.
%> @retval  V2       Scalar @f$y@f$-velocity of the current time @f$[ \#V+\#T\times 1 ]@f$ for force term e_2.
%> @retval  P        Scalar pressure of the old time step    @f$[ \#V\times 1 ]@f$ for force term e_2.


function [U, V, U2, V2, P, bU] = solveSystem(this, Uold, Vold, tau, dataN, datauD, datavD, isVisPattern)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ASSERTIONS AND ABBREVS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

assert(nargin == 8, 'HyPHM: Wrong number of arguments.')
assert(isa(this, 'Perm'))

g = this.grid;
numV = g.numV;
numE = g.numE;
numT = g.numT;
numVT = numV + numT; % abbrevation, often used

Uold = zeros(numVT, 1);
Vold = zeros(numVT, 1);
datauD = zeros(numVT, 1);
datavD = zeros(numVT, 1);
assert(length(Uold) == numVT)
assert(length(Vold) == numVT)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  markDirN: [\#V+\#E x 1 logical]

%% DEFINITION OF INDEX SETS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DIRICHLET NODES (for velocity) % for old implementation see F1
Level = this.L.getdata(this.stepper.curstep-1);

%markDirE = false(numE, 1); % [\#E x 1] % fetch edge ids acc. to Dirichlet boundary


markDirV = false(numV, 1); % [\#V x 1]

for node = 1:this.grid.numV %usual nodes in solid
    markDirV(node) = (Level(node) >= -eps);
end

if all(~markDirV) %at least one Dirichlet boundary
    markDirV(1) = true;
end
ActiveTriangles = ~markDirV(g.V0T(:, 1)) | ~markDirV(g.V0T(:, 2)) | ~markDirV(g.V0T(:, 3));
markDirE = markDirV(g.V0E(:, 1)) & markDirV(g.V0E(:, 2));

markDirN = [markDirV; ~ActiveTriangles]; % mark all nodes (offset for edge nodes)

ActiveTriangles = find(ActiveTriangles);

markNeumE = false(numE, 1); % [\#E x 1] % fetch edge ids acc. to Neumann boundary

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Pressure DOFs%%
%pressure DOF if velocity DOF or neighboring one

c = [(~markDirV(g.V0E(:, 1))) .* g.V0E(:, 2); (~markDirV(g.V0E(:, 2))) .* g.V0E(:, 1)];
c = unique(c);
if c(1) == 0
    c = c(2:end);
end

pressureDOF = ~markDirV;
pressureDOF(c(:), 1) = true;

index = find(pressureDOF); %pressure normalisation
if length(index) > 0
    pressureDOF(index(1)) = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% UNIQUENESS OF PRESSURE %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~any(markNeumE) && isempty(this.balanceP)
    %%printline((-1), 'No Neum-bdry => balanceP on!')
elseif any(markNeumE) && ~isempty(this.balanceP)
    %%printline((-1), 'balanceP on => No Neum-bdry!')
end


%%%%%%%%%%%%%%

%% ASSEMBLY %%
%%%%%%%%%%%%%%

tAssembly = tic;
%%printline(2, 'Assembly of linear terms for a Stokes problem')

% preparing data
uD = markDirN .* datauD;
vD = markDirN .* datavD; % Dirichlet data uD [#V+#T, 1], where uD has zero-entries if node is not on Gamma_D

%% Assembly of A, B, C, E
%%printline(3, 'Assembling the block matrices.')
ForceE = this.assembleSystem(ActiveTriangles); % assemble or load if already assembled

%% Assembly of the right hand-side

% for jj = 1 : 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tic
% bU and bV
%printline(3, 'Assembling bU and bV')
bU = zeros(numVT, 1); % This is build over all indices.
bV = zeros(numVT, 1); % Later, the system is solved on a sub index set only.
% force term

% for i=1:length(ActiveTriangles)
%   kT = ActiveTriangles(i);
%   nodesOnVerts = g.V0T(kT, :); % vertex numbers on T used for indexing
%   nodesOnEdges = g.E0T(kT, :); % edge numbers on T used for indexing
%   idxsU = [nodesOnVerts, numV + nodesOnEdges]; % using numV as offset
%
%   bU(idxsU) = bU(idxsU) + ForceE{kT}(:) .* Fxvec(idxsU);
%   bV(idxsU) = bV(idxsU) + ForceE{kT}(:) .* Fyvec(idxsU);
% end
bU = bU + ForceE;
% toc
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tic
% %printline(3, 'new 1')
% bbU = zeros(numVE, 1); % This is build over all indices.
% bbV = zeros(numVE, 1); % Later, the system is solved on a sub index set only.
%
% loccE = ...
%    [ 1/60,   -1/360, -1/360, -1/90,   0,     0;
%   -1/360,   1/60,  -1/360,   0,   -1/90,   0;
%   -1/360,  -1/360,  1/60,    0,     0,   -1/90;
%   -1/90,     0,      0,     4/45,  2/45,  2/45;
%   0,     -1/90,    0,     2/45,  4/45,  2/45;
%   0,       0,    -1/90,   2/45,  2/45,  4/45];
%
% % force term
% for kT = 1 : numT
%   idxsU = [g.V0T(kT, :), numV + g.E0T(kT, :)];
%   detDF = det(squeeze(g.A(kT,:,:)));
%
%   bbU(idxsU) = bbU(idxsU) + detDF * loccE * Fxvec(idxsU);
%   bbV(idxsU) = bbV(idxsU) + detDF * loccE * Fyvec(idxsU);
% end
% toc
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tic
% %printline(3, 'new 2')
%
%
% detDF = zeros(numT, 1);
% for kT = 1 : numT
%   detDF(kT) = det(squeeze(g.A(kT,:,:)));
% end
% loccE = ...
%    [ 1/60,   -1/360, -1/360, -1/90,   0,     0;
%   -1/360,   1/60,  -1/360,   0,   -1/90,   0;
%   -1/360,  -1/360,  1/60,    0,     0,   -1/90;
%   -1/90,     0,      0,     4/45,  2/45,  2/45;
%   0,     -1/90,    0,     2/45,  4/45,  2/45;
%   0,       0,    -1/90,   2/45,  2/45,  4/45];
%
% idxR = [g.V0T, g.numV + g.E0T]; % [#T x 6]
% idxC = ones(g.numT, 6);
%
% valsV = zeros(g.numT, 6);
% for kT = 1 : numT
%   valsV(kT, :) = detDF(kT) * loccE * Fyvec(idxR(kT,:));
% end
%
%
% valsU = Fxvec(idxR);
% for k = 1 : 6
%   valsU(:, k) = valsU(:, k) .* detDF;
% end
%
% for kT = 1 : g.numT
%   valsU(kT, :) = loccE*valsU(kT, :)';
% end
%
%
% bbU = sparse(idxR, idxC, valsU, numVE, 1);
% bbU = full(bbU);
% bbV = sparse(idxR, idxC, valsV, numVE, 1);
% bbV = full(bbV);
%
%
% toc
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assert(norm(bV-bbV) < 1E-17, 'bbV not correct: norm(bV-bbV) = %.3e', norm(bV-bbV))
% assert(norm(bU-bbU) < 1E-17, 'bbU not correct: norm(bU-bbU) = %.3e', norm(bU-bbU))
%
%
%
% end
% close all
% error('foo')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dirichlet term
%bU = bU - dataN*this.globA * uD;
%bV = bV - dataN*this.globA * vD;

% bP
%printline(3, 'Assembling bP')
bP = -this.globB * uD - this.globC * vD;

% cU and cV (if isEvolution)
if this.isEvolution
    %printline(3, 'Assembling cU and cV')
    cU = this.globE * Uold - this.globE * uD;
    cV = this.globE * Vold - this.globE * vD;
else
    cU = zeros(numVT, 1);
    cV = zeros(numVT, 1);
end

%% Pressure Balance Contraint               (pressure balance, row vector)
if ~isempty(this.balanceP)
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
    rhsP = sum(g.areaT) * this.balanceP;
else rhsP = 0;
end

%% Solution Vector
% write Dirichlet values to solution vector (these entries wont be touched by solver)
X = [uD; vD; zeros(numV, 1); 0]; % [U, V, P, dum] where dum is a scalar dummy for the mass balance
X = [X, X];
%printline(3,  '...done [%.3f sec]', toc(tAssembly)) % stop timer assembly linear part


% % Build large system of equations
% Z  = zeros(1, numVT);
% bigA = [ dataN*this.globA ,                      sparse(numVT,numVT),           -this.globB',              Z';
%         sparse(numVT,numVT),                     dataN*this.globA,              -this.globC',              Z';
%                this.globB,                       this.globC,                     sparse(numV,numV),         mP';
%                Z,                                Z,                              mP,                        1 ]; % additional line for pressure balance
% % clear A B C E
%
% % Build large right-hand side
% bigB1 = [  bU ;
%           bV ;
%           bP;
%           rhsP  ]; % additional entry for pressure balance

%% Solving the large system
% define subindices (solve without Dirichlet indices)
%idx = logical([~markDirN; ~markDirN; ~markDirV; mPidx]);
idx = logical([~markDirN; ~markDirN; pressureDOF; mPidx]); % mPidx = false if there's no mass balance
%printline(2, 'Number of non-zeros of system matrix');
%printline(3, '%d (%d on degrees of freedoms)', nnz(bigA), nnz(bigA(idx, idx)))

% calculate secound row of permeabilty tensor
% bigB2 = [  bV ;
%           bU;
%           bP;
%           rhsP  ];

tSolving = tic; %printline(2, 'Solving of the large system for a Stokes problem')
A = dataN * this.globA;
AS = A(~markDirN, ~markDirN);
B = [this.globB, this.globC];
B = B(pressureDOF, [~markDirN; ~markDirN]);
f = [bU; bV];
f = f([~markDirN, ~markDirN]);
fp = bP(pressureDOF);
n = sum(~markDirN) * 2;
%%%%%%%%%%%%
% Counteract scaling behavior of different matrix blocks
scale = sparse(g.V0T, 1, kron([1, 1, 1], g.areaT), g.numV, 1);
scale = 1 ./ sqrt(scale(pressureDOF));

%%%%%%%%%%%%
MY = X;
[MY(idx, 1), iter1] = Uzawa(AS, B, f, fp, zeros(sum(pressureDOF), 1), n, scale);
f = [bV; bU];
%iter1
f = f([~markDirN, ~markDirN]);
[MY(idx, 2), iter2] = Uzawa(AS, B, f, fp, zeros(sum(pressureDOF), 1), n, scale);

%%%%%%%%%%%%%%%%%%%%%%%
% tic
% options.droptol = 0.02;
% [PREC,options]=AMGfactor(bigA(idx, idx),options);
%
% options.tol=10^(-5);
% options.nrestart=1000;
% options.maxiter=100000;
% [x,options]=AMGsolver(bigA(idx, idx),PREC,options,bigB1(idx));
% [x,options]=AMGsolver(bigA(idx, idx),PREC,options,bigB1(idx));
% options.niter
% toc
%%%%%%%%%%%%%%%%%%%%%%%%

%X(idx,:) = bigA(idx, idx)\[bigB1(idx), bigB2(idx)];
%norm(MY(1 : 2*numVE,1)-X(1 : 2*numVE,1))/norm(X(1 : 2*numVE),1)
X = MY;
% extract unknowns from solution vector, set new iterates
% and compute the error
U = X(1:numVT, 1);
V = X(numVT+1:2*numVT, 1);
P = X(2*numVT+1:2*numVT+numV, 1); % there is the scalar mass balance dummy at the end of X
%printline(3,  '...done [%.3f sec]', toc(tSolving)) % stop timer solving
U2 = X(1:numVT, 2);
V2 = X(numVT+1:2*numVT, 2);

%printline(3,  '...done [%.3f sec]', toc(tSolving)) % stop timer solving
% optionally visualize sparsity pattern
if isVisPattern % visualize sparsity pattern
    figure % open new figure
    spy(bigA(idx, idx))
    title('Sparsity Pattern of System Matrix (DOFs only)')
    nz = nnz(bigA(idx, idx));
    xlabel(sprintf('non-zeros=%d (%.3f%%)', nz, nz * (100 / numel(bigA(idx, idx)))));
end

% test uniqueness
% dummytol = 1E-13;
% if abs(X(end)) > dummytol
%   error('HyPHM: The resulting system is over determined.  Remove constraints!  (...or if dummy ~ tol decrease mesh size) [dummy = %.3e, tol = %.3e]', X(end), dummytol);
% end


%%%%%%%%%%%%%%%%%%%%

%% POSTPROCESSING %%
%%%%%%%%%%%%%%%%%%%%

%% Mean pressure
% meanP = 0.0;
% for kT = 1 : numT
%   if sum(isnan(P(g.V0T(kT, :)))) == 0
%     meanP = meanP + g.areaT(kT)/3*sum(P(g.V0T(kT, :)));
%   end
% end
index = sum(isnan(P(g.V0T(:, :))), 2) == 0;
meanP = sum(g.areaT(index)./3.*sum(P(g.V0T(index, :)), 2));
meanP = meanP / sum(g.areaT);
%printline(3, 'Mean scalar solution is %e', meanP)

return

end


%%%%%%%%%%%%%

%% REMARKS %%
%%%%%%%%%%%%%
%
% The speedup for allocating the sparse matrices is only about 1% => no
% allocation
