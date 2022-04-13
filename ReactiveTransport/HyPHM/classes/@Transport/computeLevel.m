%> @file Transport/computeLevel.m Compute the next time level for the unknowns @f$\vec{q}@f$ and @f$u@f$.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @brief Compute the next time level for the unknowns @f$\vec{q}@f$ and @f$u@f$.
%>
%> @param this       [ Transport ]  contains data required by solver
%> @param varargin   [ list ]  <code>'s'</code> or <code>'silent'</code>: silent mode, i.e., no output
%> @param varargin   [ list ]  <code>'spy'</code>: visualize the sparsity pattern of the system matrix
%> @param varargin   [ list ]  <code>'BDF2'</code> Use BDF2 for time discretization <b>(only equal time steps!)</b>
%>
%> @sa HyPHM, Grid, Stepper, AbstractProblem
%>


function [bigA, bigB] = computeLevel(this, varargin)

if nargin > 1 % if varargin is non-empty
    if ~ischar(varargin{1}) % if old call solve(this, Variable, Variable) is used
        msg = 'HyPHM: SYNTAX CHANGE: Please include Your Unknowns for mass flux and concentration in the class Transport by setting Transport.Q = ..., Transport.U = ... and call Transport.computeLevel() without unknowns.';
        error(msg);
    end
end

%% Options.
if any(ismember(varargin, 's')) || any(ismember(varargin, 'silent'))
    isSlt = true;
else
    isSlt = false;
end
if any(ismember(varargin, 'spy'))
    isVisPattern = true;
else
    isVisPattern = false;
end
if any(ismember(varargin, 'BDF2'))
    BDFq = 2; % BDF2 for time discretization.
else
    BDFq = 1; % BDF1 = impl. Euler for time discretization.
end

st = this.stepper;
numT = this.grid.numT;

%% Checking Coefficient Functions.
initVoidCoeffs(this, isSlt) % initialize non-defined coeffs by zeros

%% Preparing Input.
Upasts{1} = this.U.getdata(st.curstep-1);
if BDFq > 1, Upasts{2} = this.U.getdata(st.curstep-2); end
curTau = st.curtau;

%% Preparing raw Data.
dataA = this.A.getdata(st.curstep);
dataAold = this.A.getdata(st.curstep-1);
if BDFq > 1, dataAvold = this.A.getdata(st.curstep-2);
else dataAvold = NaN;
end
dataB = this.B.getdata(st.curstep);
dataBold = this.B.getdata(st.curstep-1);
dataC = this.C.getdata(st.curstep);
dataD = this.D.getdata(st.curstep);

%% Check weather D is stationary (system matrix for D can be reused in this case).
if isequal(dataD, this.D.getdata(st.curstep - 1))
    isDstationary = true;
else
    isDstationary = false;
end


if strcmp(this.D.type, 'CCCC') % if D is in CCCC, not P0P0P0P0
    dataD_P0P0P0P0 = zeros(numT, 4);
    row = reshape(dataD, 1, 4);
    for kT = 1:numT
        dataD_P0P0P0P0(kT, :) = row;
    end
    dataD = dataD_P0P0P0P0;
else
    assert(strcmp(this.D.type, 'P0P0P0P0'), 'HyPHM: Coefficient D has to have type CCCC or P0P0P0P0 but is %s.', this.D.type)
end
assert(isequal(size(dataD), [numT, 4]), 'HyPHM kernel: Some strange things happen here.')


dataE = this.E.getdata(st.curstep);
dataF = this.F.getdata(st.curstep);
datauD = this.uD.getdata(st.curstep);
RobinData = this.uR.getdata(st.curstep);
datagF = this.gF.getdata(st.curstep);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check positive Definiteness of D.
msg = sprintf('HyPHM: Coefficient D is not positive definite in problem ''%s''!', this.name);
for kT = 1:numT
    DT = dataD(kT, :);
    if ~(DT(1) * DT(4) - DT(2) * DT(3) > 0 && DT(1) > 0)
        error([msg, sprintf(' D at triangle %d was [%.3e, %.3e; %.3e, %.3e]', kT, DT(1), DT(2), DT(3), DT(4))])
    end
end

%% Definition of Index Sets.
[markDirE,markNeumE,markFluxE,markRobinE,markFreeE,idxDirE,idxNeumE,idxFluxE,idxRobinE,idxFreeE] = ...
  tp_indexsets(this);

%% Assembly.
[B, C, D, E, bQ, bU, Y] = ...
  assembleSystem(this, Upasts, BDFq, curTau, ...
  dataA, dataAold, dataAvold, dataB, dataBold, dataC, dataD, RobinData, isDstationary, dataE, dataF, ...
  datauD, datagF, markDirE, markNeumE, markRobinE,markFluxE, idxNeumE, idxRobinE, isSlt); 

%% Calling the linear Solver.
[Qnew, Unew, bigA, bigB] = ...
    solveSystem(this, curTau, B, C, D, E, bQ, bU, markNeumE, markFreeE, dataC, dataE, Y, isSlt, isVisPattern);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Reprocessing Output
this.Q.setdata(st.curstep, Qnew);
this.U.setdata(st.curstep, Unew);

end
