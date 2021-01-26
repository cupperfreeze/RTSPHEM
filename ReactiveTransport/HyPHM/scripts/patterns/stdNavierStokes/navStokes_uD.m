% Copyright 2009, 2010 F. Frank, Chair of Applied Mathematics I,
%        Department of Mathematics, University of Erlangen-Nuremberg,
%        G E R M A N Y

function ret = navStokes_uD(~, x)

assert(nargin == 2)
assert(isequal(size(x), [2, 1]), 'x has wrong dimensions')

ret = [0; 0];

end
