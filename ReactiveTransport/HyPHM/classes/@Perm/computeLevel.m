%> @file Perm/computeLevel.m Compute the next time level for the unknowns @f$\vec{u}@f$ and @f$p@f$.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @brief Compute the unknowns @f$\vec{u}@f$ and @f$p@f$..
%> This routine is a specialized version of 'Stokes' to effectively solve
%> permeability cell problems. For any other purpose, use 'Stokes'

%>
%> Further, this wrapper initializes the void coefficient functions by zero or one (viscosity).
%>
%> @param this Instance of Perm
%> @param varargin   [ list ]  's' or 'silent': silent mode, ie no output
%>
%> @sa HyPHM, Stokes, Transport
%> Returns integrals of basis functions in bU [numV+numT x 1] or [numV+numE x 1]


function [bU] = computeLevel(this, varargin)

if nargin > 1 % if varargin is non-empty

    if ~ischar(varargin{1}) % if old call solve(this, Variable, Variable) is used
        msg = 'HyPHM: SYNTAX CHANGE: Please include Your Unknowns for velocity and pressure in the class Stokes by setting Stokes.U = ..., Stokes.P = ... and call Stokes.computeLevel() without unknowns.';
        error(msg);
    end
end

%% Options
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

%% Check Unknowns and Coefficient Functions
initVoidCoeffs(this) % initialize non-defined coeffs by zeros

%% Preparing Input
st = this.stepper;
Uvecold = this.U.getdata(st.curstep-1);
Uolddata = Uvecold(:, 1);
Volddata = Uvecold(:, 2);

%dataF  = this.F.getdata(st.curstep);
dataN = this.N.getdata(st.curstep);
datauDvec = this.uD.getdata(st.curstep);
datauD = datauDvec(:, 1);
datavD = datauDvec(:, 2);

%% Calling Solver depending on chosen discretization
if this.Bubble
    [Unew, Vnew, U2new, V2new, Pnew, bU] = solveSystemBubble(this, Uolddata, Volddata, st.curtau, dataN, datauD, datavD, isVisPattern);
else
    [Unew, Vnew, U2new, V2new, Pnew, bU] = solveSystem(this, Uolddata, Volddata, st.curtau, dataN, datauD, datavD, isVisPattern);
end

%% Meging Output
this.U.setdata(st.curstep, [Unew, Vnew]);
this.U2.setdata(st.curstep, [U2new, V2new]);
this.P.setdata(st.curstep, Pnew);

end
