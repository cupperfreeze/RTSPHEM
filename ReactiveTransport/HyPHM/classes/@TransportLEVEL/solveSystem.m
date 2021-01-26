%> @file TransportLEVEL/solveSystem.m Solves the linear system of the Transport problem

function [Qnew, Unew, bigA, bigB] = solveSystem(d, curTau, B, C, D, E, bQ, bU, markNeumE, markFreeE, dataC, dataE, Y, isSlt, isVisPattern, markDirT)

g = d.grid;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Concentration Balance Contraint %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% defining index for mass balance mUidx
if ~isempty(d.balanceU);
    mUidx = true;
else mUidx = false;
end % mUidx: mean U index {0,1}

% assembly of mU (column vector for mass balance)
if mUidx;
    mU = g.areaT;
else mU = zeros(g.numT, 1);
end

% building right-hand side rhsU
if mUidx;
    rhsU = sum(g.areaT) * d.balanceU;
else rhsU = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Large System of Equations %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Z = zeros(g.numE, 1);
bigA = [B, C + D'; ...
    curTau * D, E; ...
    Z', mU']; clear B C D E mU % remove large matrices from memory

bigB = [bQ; bU; rhsU]; clear bQ bU % remove large vectors from memory
% define subindices (solve without Dirichlet indices)
if sum(g.idE) == 0;
    noDirT = find(~markDirT);
    markDirT(noDirT(1)) = true;
end
idx1 = logical([markFreeE; ~markDirT; mUidx]); % mUidx = false if there's no mass balance
idx2 = logical([markFreeE; ~markDirT]);
X = [Y; zeros(g.numT, 1)];
printline(~isSlt*2, 'Number of non-zeros of system matrix');
printline(~isSlt*3, '%d (%d on degrees of freedoms)', nnz(bigA), nnz(bigA(idx1, idx2)))
% optionally visualize sparsity pattern
if isVisPattern % visualize sparsity pattern
    figure % open new figure
    spy(bigA(idx1, idx2))
    title('Sparsity Pattern of System Matrix (DOFs only)')
    nz = nnz(bigA(idx1, idx2));
    xlabel(sprintf('non-zeros=%d (%.3f%%)', nz, nz * (100 / numel(bigA(idx1, idx2)))));
end
% start timer
printline(~isSlt*2, 'Solving of the large system for the problem ''%s''', d.name)
tSolving = tic;
% call linear solver
try % this may run out of memory

    X(idx2) = bigA(idx1, idx2) \ bigB(idx1);

    %   tic
    %   M = bigA(idx1, idx2);
    %   b =  bigB(idx1);
    %
    %   options.droptol=10^(-3);
    %   [PREC,options]=AMGfactor(M,[]);
    %   options.nrestart = 100;
    %   options.maxit=10000;
    %   options.restol=1e-10;
    %   [sol, options] = AMGsolver(M, PREC, options, b);
    %   PREC = AMGdelete(PREC);
    %   clear options AMGfactor AMGsolver
    %   X(idx2) = sol;
    %   toc
    %   norm(X-Y)
    %   norm(X)
    %  X=Y;
catch exception
    printline(1, 'Direct solver run out of memory, using iterative solver lsqr instead (time consuming).')
    lsqr_tol = 1E-6;
    lsqr_maxiter = 1E6;
    printline(3, 'Trying to solve with tol = %.3e and maxiter = %d', lsqr_tol, lsqr_maxiter)
    X(idx2) = lsqr(bigA(idx1, idx2), bigB(idx1), lsqr_tol, lsqr_maxiter);
    printline(3, 'lsqr did succeed!')
end

%% Remark 4: Building a quadratic System (does not work properly, see below)

printline(~isSlt*3, '                      ...done [%.3f sec]', toc(tSolving)) % stop timer solving

%%%%%%%%%%%%%%%%%%%%

%% Postprocessing %%
%%%%%%%%%%%%%%%%%%%%

%% Extraction of representation vectors
Qnew = X(1:g.numE);
Unew = X(g.numE+1:g.numT+g.numE);

%% Additional flux terms for Neumann edges (remark 5)
Qnew(markNeumE) = Qnew(markNeumE) + Unew(g.T0E(markNeumE, 1)) .* dataC(markNeumE) + dataE(markNeumE);

%% Mean concentration
meanU = dot(Unew, g.areaT) / sum(g.areaT);
printline(~isSlt*3, 'Mean scalar solution is %f', meanU)

end

%%%%%%%%%%%%%%%

%% Footnotes %%
%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Remark 3: Solving via ILUpack
% case {'ilupack'}
%   path(path, './opt/ilupack/matlab')
%   options.isdefinite = 0;
%   options.amg = 'ilu'; % 'amli'  'mg'
%   [PREC, options]        = AMGfactor(A(idx, idx));
%   [x(idx), options]   = AMGsolver(A(idx, idx), PREC, options, b(idx), x(idx, T-1)); %#ok<NASGU>
%   PREC                   = AMGdelete(PREC);                      %#ok<NASGU> % free memory

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Remark 4: Building a quadratic System (does not work properly)
% Z = zeros(numE, 1);
% bigA = [ B,   C,   Z;
%          D,   E,  mU;
%          Z', mU',  1];      clear B C D E mU
% bigB = [bQ; bU;  rhsU];     clear bQ bU
% % define subindices (solve without flux and Neumann indices)
% idx = logical([markFreeE; ones(numT, 1); mUidx]); % mUidx = false if there's no mass balance
% tSolving = tic; printline(~isSlt*2, 'Solving of the large system for the Transport problem')
% X = [y; zeros(numT, 1); 0];
% X(idx) = bigA(idx, idx) \ bigB(idx);
% % test uniqueness
% if abs(X(end)) > 1E-15
%   error('HyPHM: The resulting system is over determined.  Remove constraints! [dummy = %.3e]', X(end));
% end