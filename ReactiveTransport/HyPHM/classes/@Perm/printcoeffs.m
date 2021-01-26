% Copyright 2009, 2010 F. Frank, Chair of Applied Mathematics I,
%        Department of Mathematics, University of Erlangen-Nuremberg,
%        G E R M A N Y

function printcoeffs(obj)

printline(2, 'Coefficient functions for Stokes');

printline(3, 'Kinematik viscosity N:  %s', x2str(obj.N));
printline(3, 'Source/sink         F:  %s', x2str(obj.F));
printline(3, 'Dirichlet          uD:  %s', x2str(obj.uD));

end
