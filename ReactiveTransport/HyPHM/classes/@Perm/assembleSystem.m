%> @file Perm/assembleSystem.m Assemble the linear system for the Stokes problem.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ForceE = assembleSystem(this, ActiveTriangles)

% if A, B, C, E are already assembled

% if ~isempty(this.globA) && ~isempty(this.globB) && ~isempty(this.globC) && ~isempty(this.globE)
%
%   return % do nothing
%   % if none of the matrices is assembled yet, compute them
% elseif isempty(this.globA) && isempty(this.globB) && isempty(this.globC) && isempty(this.globE)
% restrict assembly to traingles in fluid domain via 'ActiveTriangles'

numV = this.grid.numV;
numE = this.grid.numE;
numT = length(ActiveTriangles); %this.grid.numT;
numVE = numV + numE;
numVT = this.grid.numT + numV;
g = this.grid;


%   Arow = cell(numT, 1);  Acol = cell(numT, 1);  Aval = cell(numT, 1);
%   Brow = cell(numT, 1);  Bcol = cell(numT, 1);  Bval = cell(numT, 1);
%   Crow = cell(numT, 1);  Ccol = cell(numT, 1);  Cval = cell(numT, 1);
%   Erow = cell(numT, 1);  Ecol = cell(numT, 1);  Eval = cell(numT, 1);
%
%   ForceE = cell(numT, 1);
%   %Note.  Check for parallel toolbox via "if PCTInstalledAndLicensed".
%   for i = 1 : numT
%     kT = ActiveTriangles(i);
%     nodesOnVerts = g.V0T(kT, :); % vertex numbers on T used for indexing
%     nodesOnEdges = g.E0T(kT, :); % edge numbers on T used for indexing
%     idxsU = [nodesOnVerts, numV + nodesOnEdges]; % using numV as offset
%     idxsP = nodesOnVerts;
%     DF = reshape(g.A(kT,:,:),2,2); % the jacobian of the affine map F
%
%      Arow{kT} = repmat(idxsU', 1, 6);  % all cells are 6 x 6
%      Acol{kT} = repmat(idxsU,  6, 1);
%     Aval{kT} = this.locA(DF);
%
%      Brow{kT} = repmat(idxsP', 1, 6);  % all cells are 3 x 6
%      Bcol{kT} = repmat(idxsU,  3, 1);
%     Bval{kT} = this.locB(DF);
%
%      Crow{kT} = repmat(idxsP', 1, 6);  % all cells are 3 x 6
%      Ccol{kT} = repmat(idxsU,  3, 1);
%    Cval{kT} = this.locC(DF);
%
%     ForceE{kT} = this.locE(DF);
%
%    end
%%%%%%%%%%%%%%%%%%%%                                                           vectorized assembly by SG
if ~this.Bubble
    nodesOnVerts = g.V0T(ActiveTriangles, :); % vertex numbers on T used for indexing
    nodesOnEdges = g.E0T(ActiveTriangles, :); % edge numbers on T used for indexing
    idxsU = [nodesOnVerts, numV + nodesOnEdges]; % using numV as offset
    idxsP = nodesOnVerts;
    DF = g.A(ActiveTriangles, :, :); % the jacobian of the affine map F

    Arow = reshape(repmat(idxsU', [1, 1, 6]), [numT * 6, 6]);
    Acol = repmat(idxsU, [1, 1, 6]);
    Acol = permute(Acol, [3, 1, 2]);
    Acol = reshape(Acol, numT*6, 6);
    Aval = this.locA(DF, numT);
    Aval = permute(Aval, [1, 3, 2]);
    Aval = reshape(Aval, numT*6, 6);

    Brow = reshape(repmat(idxsP', [1, 1, 6]), [numT * 6, 3]);
    Bcol = repmat(idxsU, [1, 1, 3]);
    Bcol = permute(Bcol, [3, 1, 2]);
    Bcol = reshape(Bcol, numT*6, 3);
    Bval = this.locB(DF, numT);
    Bval = permute(Bval, [1, 3, 2]);
    Bval = reshape(Bval, numT*6, 3);

    Cval = this.locC(DF, numT);
    Cval = permute(Cval, [1, 3, 2]);
    Cval = reshape(Cval, numT*6, 3);

    Frow = idxsU;
    Fcol = 1;
    Fval = this.locE(DF, numT);
    Fval = permute(Fval, [1, 3, 2]);
    Fval = reshape(Fval, numT*6, 1);
    ForceE = sparse(Frow, Fcol, Fval, numVE, 1);
    try
        A = fsparse(Arow, Acol, Aval, [numVE, numVE], [], 1);
        B = fsparse(Brow, Bcol, Bval, [numV, numVE], [], 1);
        C = fsparse(Brow, Bcol, Cval, [numV, numVE], [], 1);
    catch
        A = sparse(Arow, Acol, Aval, numVE, numVE);
        B = sparse(Brow, Bcol, Bval, numV, numVE);
        C = sparse(Brow, Bcol, Cval, numV, numVE);

    end

else
    nodesOnVerts = g.V0T(ActiveTriangles, :); % vertex numbers on T used for indexing
    nodesBary = ActiveTriangles; % edge numbers on T used for indexing
    idxsU = [nodesOnVerts, numV + nodesBary]; % using numV as offset
    idxsP = nodesOnVerts;
    DF = g.A(ActiveTriangles, :, :); % the jacobian of the affine map F
    Arow = reshape(repmat(idxsU', [1, 1, 4]), [numT * 4, 4]);
    Acol = repmat(idxsU, [1, 1, 4]);
    Acol = permute(Acol, [3, 1, 2]);
    Acol = reshape(Acol, numT*4, 4);
    Aval = this.locABubble(DF, numT);
    Aval = permute(Aval, [1, 3, 2]);
    Aval = reshape(Aval, numT*4, 4);

    Brow = reshape(repmat(idxsP', [1, 1, 4]), [numT * 4, 3]);
    Bcol = repmat(idxsU, [1, 1, 3]);
    Bcol = permute(Bcol, [3, 1, 2]);
    Bcol = reshape(Bcol, numT*4, 3);
    Bval = this.locBBubble(DF, numT);
    Bval = permute(Bval, [1, 3, 2]);
    Bval = reshape(Bval, numT*4, 3);

    Cval = this.locCBubble(DF, numT);
    Cval = permute(Cval, [1, 3, 2]);
    Cval = reshape(Cval, numT*4, 3);

    Frow = idxsU;
    Fcol = 1;
    Fval = this.locEBubble(DF, numT);
    Fval = permute(Fval, [1, 3, 2]);
    Fval = reshape(Fval, numT*4, 1);
    ForceE = sparse(Frow, Fcol, Fval, numVT, 1);
    try
        A = fsparse(Arow, Acol, Aval, [numVT, numVT], [], 1);
        B = fsparse(Brow, Bcol, Bval, [numV, numVT], [], 1);
        C = fsparse(Brow, Bcol, Cval, [numV, numVT], [], 1);
    catch
        A = sparse(Arow, Acol, Aval, numVT, numVT);
        B = sparse(Brow, Bcol, Bval, numV, numVT);
        C = sparse(Brow, Bcol, Cval, numV, numVT);

    end


end


%%%%%%%%%%%%%%%%%%%%%

%A = sparse(cell2mat(Arow), cell2mat(Acol), cell2mat(Aval), numVE, numVE);
%  B = sparse(cell2mat(Brow), cell2mat(Bcol), cell2mat(Bval), numV,  numVE);
% C = sparse(cell2mat(Crow), cell2mat(Ccol), cell2mat(Cval), numV,  numVE);
if this.isEvolution % otherwise equal to zero
    E = sparse(cell2mat(Erow), cell2mat(Ecol), cell2mat(Eval), numVE, numVE);
else
    E = sparse(numVE, numVE);
end
% clear Arow Acol Aval Brow Bcol Bval Crow Ccol Cval Erow Ecol Eval Fval
% See Footnote 1 below for the naive assemblies.

this.globA = A;
this.globB = B;
this.globC = C;
this.globE = E;


% else % if some thing are strange
%   error('HyPHM: Something strange happens in the kernel.')
% end

end

%% Footnote 1 (Naive assemblies).
%   %% NAIVE SINGLE CORE VERSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   A = sparse(numVE, numVE); % quadratic sparse matrix which is build over all indices. Later, the system is solved on a sub index set only.
%   B = sparse(numV,  numVE);
%   C = sparse(numV,  numVE);
%   E = sparse(numVE, numVE); % only assembled if non-stationary, otherwise equal to zero
%     for kT = 1 : numT
%       nodesOnVerts = g.V0T(kT, :); % vertex numbers on T used for indexing
%       nodesOnEdges = g.E0T(kT, :); % edge numbers on T used for indexing
%       idxsU = [nodesOnVerts, numV + nodesOnEdges]; % using numV as offset
%       idxsP = nodesOnVerts;
%
%       DF = squeeze(g.A(kT,:,:)); % the jacobian of the affine map F
%
%       A(idxsU, idxsU) = A(idxsU, idxsU) + this.locA(DF);
%       B(idxsP, idxsU) = B(idxsP, idxsU) + this.locB(DF);
%       C(idxsP, idxsU) = C(idxsP, idxsU) + this.locC(DF);
%
%       if this.isEvolution % otherwise equal to zero
%         E(idxsU, idxsU) = E(idxsU, idxsU) + this.locE(DF);
%       end
%     end
%     %%% NAIVE SINGLE CORE VERSION END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  %% NAIVE PARALLEL VERSION (does not scale properly) %%%%%%%%%%%%%%%%%%%%
%   A = sparse(numVE, numVE); % quadratic sparse matrix which is build over all indices. Later, the system is solved on a sub index set only.
%   B = sparse(numV,  numVE);
%   C = sparse(numV,  numVE);
%   E = sparse(numVE, numVE); % only assembled if non-stationary, otherwise equal to zero
%     spmd
%       A_ = sparse(numVE, numVE); % quadratic sparse matrix which is build over all indices. Later, the system is solved on a sub index set only.
%       B_ = sparse(numV, numV + numE);
%       C_ = sparse(numV, numV + numE);
%       E_ = sparse(numVE, numVE);
%       for kT = floor((numT/numlabs)*(labindex-1)+1) : floor((numT/numlabs)*(labindex))
%         nodesOnVerts = g.V0T(kT, :); % vertex numbers on T used for indexing
%         nodesOnEdges = g.E0T(kT, :); % edge numbers on T used for indexing
%         idxsU = [nodesOnVerts, numV + nodesOnEdges]; % using numV as offset
%         idxsP = nodesOnVerts;
%
%         DF = squeeze(g.A(kT,:,:)); % the jacobian of the affine map F
%
%         A_(idxsU, idxsU) = A_(idxsU, idxsU) + this.locA(DF);
%         B_(idxsP, idxsU) = B_(idxsP, idxsU) + this.locB(DF);
%         C_(idxsP, idxsU) = C_(idxsP, idxsU) + this.locC(DF);
%
%         if this.isEvolution % otherwise equal to zero
%           E_(idxsU, idxsU) = E_(idxsU, idxsU) + this.locE(DF);
%         end
%       end
%     end % spmd
%     for kP = 1 : numVTgth(A_)
%       A = A + A_{kP}; B = B + B_{kP}; C = C + C_{kP}; E = E + E_{kP};
%     end % collecting subdomain entries
%     clear A_ B_ C_ E_ kT
%     %%% NAIVE PARALLEL VERSION END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
