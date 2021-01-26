%> @file Perm/initVoidCoeffs.m Initialize non-defined coefficients.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @param this Instance of Perm
%> @param isSlt Flag for silent mode.

function initVoidCoeffs(this, isSlt)
numE = this.grid.numE;
numV = this.grid.numV;
if this.Bubble
    numT = this.grid.numT;
    if isempty(this.U)
        %  printline(~isSlt*3, 'Unknown U (velocity) was not defined, initializing it by zero for the initial step.')
        this.U = Variable(this.grid, this.stepper, 'U', 'P1BubbleP1Bubble');
        this.U.setdata(0, sparse(numT + numV, 2));
    end
    if isempty(this.U2)
        %  printline(~isSlt*3, 'Unknown U2 (velocity) was not defined, initializing it by zero for the initial step.')
        this.U2 = Variable(this.grid, this.stepper, 'U', 'P1BubbleP1Bubble');
        this.U2.setdata(0, sparse(numT + numV, 2));
    end
else
    if isempty(this.U)
        %  printline(~isSlt*3, 'Unknown U (velocity) was not defined, initializing it by zero for the initial step.')
        this.U = Variable(this.grid, this.stepper, 'U', 'P2P2');
        this.U.setdata(0, sparse(numE + numV, 2));
    end
    if isempty(this.U2)
        %  printline(~isSlt*3, 'Unknown U2 (velocity) was not defined, initializing it by zero for the initial step.')
        this.U2 = Variable(this.grid, this.stepper, 'U', 'P2P2');
        this.U2.setdata(0, sparse(numE + numV, 2));
    end
end
if isempty(this.P)
    %  printline(~isSlt*3, 'Unknown P (pressure) was not defined, initializing it by zero for the initial step.')
    this.P = Variable(this.grid, this.stepper, 'P', 'P1');
    this.P.setdata(0, sparse(numV, 1));
end
if isempty(this.L)
    %  printline(~isSlt*3, 'Unknown L (LevelSet) was not defined, initializing it by zero for the initial step.')
    this.L = Variable(this.grid, this.stepper, 'L', 'P1');
    this.L.setdata(0, sparse(numV, 1));
end
if isempty(this.N)
    %  printline(~isSlt*3, 'Coefficient N (viscosity) was not defined, initializing it by one for all steps.')
    this.N = Variable(this.grid, this.stepper, 'N', 'C');
    this.N.setdata(1.0);
end
%   if isempty(this.F)
%   %  printline(~isSlt*3, 'Coefficient F (force) was not defined, initializing it by zero for all steps.')
%     this.F  = Variable(this.grid, this.stepper, 'F',  'P2P2'); this.F.setdata([ones(numE+numV, 1), zeros(numE+numV,1)]);
%   end
if isempty(this.uD)
    %  printline(~isSlt*3, 'Coefficient uD (Dirichlet data) was not defined, initializing it by zero/no slip condition for all steps.')
    this.uD = Variable(this.grid, this.stepper, 'uD', 'P2P2');
    this.uD.setdata(sparse(numE + numV, 2));
end
end
