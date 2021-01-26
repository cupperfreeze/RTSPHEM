%> @file StokesLEVEL/computeLevel.m Compute the next time level for the unknowns @f$\vec{u}@f$ and @f$p@f$.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @brief Compute the next time level for the unknowns @f$\vec{u}@f$ and @f$p@f$..
%>
%> Solve the defined problem for @f$t_\mathrm{old}\rightarrow t_\mathrm{new}@f$,
%> where @f$t_\mathrm{new} = @f$Stepper.curtime.  First argument has to be
%> of type <tt>'P2P2'</tt> and second of <tt>'P1'</tt>.
%>
%> Further, this wrapper initializes the void coefficient functions by zero or one (viscosity).
%>
%> @param this Instance of Stokes
%> @param varargin   [ list ]  's' or 'silent': silent mode, ie no output
%>
%> @sa HyPHM, Stokes, Transport
%> @todo implement silent function for Stokes.
%> @todo implement info struct for Stokes.

function [bU] = computeLevel(this, varargin)

if nargin > 1 % if varargin is non-empty

    if ~ischar(varargin{1}) % if old call solve(this, Variable, Variable) is used
        msg = 'HyPHM: SYNTAX CHANGE: Please include Your Unknowns for velocity and pressure in the class StokesLEVEL by setting StokesLEVEL.U = ..., StokesLEVEL.P = ... and call StokesLEVEL.computeLevel() without unknowns.';
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

dataF = this.F.getdata(st.curstep);
dataN = this.N.getdata(st.curstep);
datauDvec = this.uD.getdata(st.curstep);
datauD = datauDvec(:, 1);
datavD = datauDvec(:, 2);

%% Calling Solver
[Unew, Vnew, Pnew] = solveSystem(this, Uolddata, Volddata, st.curtau, dataF, dataN, datauD, datavD, isVisPattern);

%% Meging Output
this.U.setdata(st.curstep, [Unew, Vnew]);
this.P.setdata(st.curstep, Pnew);

end
