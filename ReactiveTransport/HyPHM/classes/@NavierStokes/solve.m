%> @file NavierStokes/solve.m This is the user-friendly wrapper for the (private, hidden) NavierStokes.nssolve.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @brief This is the user-friendly wrapper for the (private, hidden) NavierStokes.stokessolve.
%>
%> Solve the defined problem for @f$t_\mathrm{old}\rightarrow t_\mathrm{new}@f$,
%> where @f$t_\mathrm{new} = @f$Stepper.curtime.  First argument has to be
%> of type <tt>'P2P2'</tt> and second of <tt>'P1'</tt>.
%>
%> Further, this wrapper initializes the void coefficient functions by zero or one (viscosity).
%>
%> @param this Instance of NavierStokes
%> @param U    Flux solution in P2P2 basis [ Variable ]
%> @param P    Scalar solution in P1 basis [ Variable ]
%> @param varargin   [ list ]  's' or 'silent': silent mode, ie no output
%>
%> @sa HyPHM, NavierStokes, Transport

function [info] = solve(this, varargin)

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

%% Checking Coefficient Functions
initVoidCoeffs(this) % initialize non-defined coeffs by zeros

%% Preparing Input
st = this.stepper;
Uvecold = this.U.getdata(st.curstep-1);
Uold = Uvecold(:, 1);
Vold = Uvecold(:, 2);

dataF = this.F.getdata(st.curstep);
dataN = this.N.getdata(st.curstep);
datauDvec = this.uD.getdata(st.curstep);
datauD = datauDvec(:, 1);
datavD = datauDvec(:, 2);

%% Calling Solver
[Unew, Vnew, Pnew, info] = nssolve(this, Uold, Vold, st.curtau, dataF, dataN, datauD, datavD, isSlt);

%% Meging Output
this.U.setdata(st.curstep, [Unew, Vnew]);
this.P.setdata(st.curstep, Pnew);

end


%> @brief Initialize non-defined coefficients by zero or one (viscosity).
%> @param this Instance of NavierStokes
function initVoidCoeffs(this)
numE = this.grid.numE;
numV = this.grid.numV;

if isempty(this.U)
    printline(3, 'Unknown U (velocity) was not defined, initializing it by zero for the initial step.')
    this.U = Variable(this.grid, this.stepper, 'U', 'P2P2');
    this.U.setdata(0, zeros(numE + numV, 2));
end
if isempty(this.P)
    printline(3, 'Unknown P (pressure) was not defined, initializing it by zero for the initial step.')
    this.P = Variable(this.grid, this.stepper, 'P', 'P1');
    this.P.setdata(0, zeros(numV, 1));
end
if isempty(this.N)
    printline(3, 'Coefficient N (viscosity) was not defined, initializing it by one.')
    this.N = Variable(this.grid, this.stepper, 'N', 'C');
    this.N.setdata(1.0);
end
if isempty(this.F)
    printline(3, 'Coefficient F (force) was not defined, initializing it by zero.')
    this.F = Variable(this.grid, this.stepper, 'F', 'P2P2');
    this.F.setdata(zeros(numE + numV, 2));
end
if isempty(this.uD)
    printline(3, 'Coefficient uD (Dirichlet data) was not defined, initializing it by zero/no slip condition.')
    this.uD = Variable(this.grid, this.stepper, 'uD', 'P2P2');
    this.uD.setdata(zeros(numE + numV, 2));
end
end
