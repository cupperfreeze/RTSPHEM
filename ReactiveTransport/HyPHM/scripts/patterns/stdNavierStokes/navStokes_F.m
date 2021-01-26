% Copyright 2009, 2010 F. Frank, Chair of Applied Mathematics I,
%        Department of Mathematics, University of Erlangen-Nuremberg,
%        G E R M A N Y

function ret = navStokes_F(~, x)

assert(isequal(size(x), [2, 1]), 'x has wrong dimensions')

if norm(x) <= 0.2
    ret = [-1; -1];
else
    ret = [0; 0];
end

end
