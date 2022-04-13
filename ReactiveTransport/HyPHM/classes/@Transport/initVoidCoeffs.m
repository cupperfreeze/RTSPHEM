%> @file initVoidCoeffs.m Initialize non-defined coefficients by zero.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @param this Instance of Transport.
%> @param isSlt Flag for silent mode.

function initVoidCoeffs(this, isSlt)
numE = this.grid.numE;
numT = this.grid.numT;

if isempty(this.Q)
    printline(~isSlt*3, 'Unknown Q was not defined, initializing it by zero for the initial step.')
    this.Q = Variable(this.grid, this.stepper, 'Q', 'RT0');
    this.Q.setdata(0, zeros(numE, 1));
end
if isempty(this.U)
    printline(~isSlt*3, 'Unknown U was not defined, initializing it by zero for the initial step.')
    this.U = Variable(this.grid, this.stepper, 'U', 'P0');
    this.U.setdata(0, zeros(numT, 1));
end
if isempty(this.A)
    printline(~isSlt*3, 'Coefficient A was not defined, initializing it by zero for all steps.')
    this.A = Variable(this.grid, this.stepper, 'A', 'P0');
    this.A.setdata(zeros(numT, 1));
end
if isempty(this.B)
    printline(~isSlt*3, 'Coefficient B was not defined, initializing it by zero for all steps.')
    this.B = Variable(this.grid, this.stepper, 'B', 'P0');
    this.B.setdata(zeros(numT, 1));
end
if isempty(this.C)
    printline(~isSlt*3, 'Coefficient C was not defined, initializing it by zero for all steps.')
    this.C = Variable(this.grid, this.stepper, 'C', 'RT0');
    this.C.setdata(zeros(numE, 1));
end
if isempty(this.D)
    printline(~isSlt*3, 'Coefficient D was not defined, initializing it by the unit matrix for all steps.')
    this.D = Variable(this.grid, this.stepper, 'D', 'CCCC');
    this.D.setdata(eye(2));
end
if isempty(this.E)
    printline(~isSlt*3, 'Coefficient E was not defined, initializing it by zero for all steps.')
    this.E = Variable(this.grid, this.stepper, 'E', 'RT0');
    this.E.setdata(zeros(numE, 1));
end
if isempty(this.F)
    printline(~isSlt*3, 'Coefficient F was not defined, initializing it by zero for all steps.')
    this.F = Variable(this.grid, this.stepper, 'F', 'P0');
    this.F.setdata(zeros(numT, 1));
end
if isempty(this.uD)
    printline(~isSlt*3, 'Dirichlet boundary data uD was not defined, initializing it by NaN''s for all steps (cautionary).')
    this.uD = Variable(this.grid, this.stepper, 'uD', 'P0E');
    this.uD.setdata(NaN(numE, 1));
end
if isempty(this.gF)
    printline(~isSlt*3, 'Flux boundary data gF was not defined, initializing it by NaN''s for all steps (cautionary).')
    this.gF = Variable(this.grid, this.stepper, 'gF', 'P0E');
    this.gF.setdata(NaN(numE, 1));
end
if isempty(this.uR)
    printline(~isSlt*3, 'Robin ambient data uR was not defined, initializing it by NaN (cautionary).')
    this.uR  = Variable(this.grid, this.stepper, 'uR',  'P0E'); this.uR.setdata(NaN(numE, 1));
end
if isempty(this.rob)
    printline(~isSlt*3, 'Robin constant rob was not defined, initializing it by NaN (cautionary).')
    this.rob  = Variable(this.grid, this.stepper, 'rob',  'C'); this.rob.setdata(NaN);
end 

end
