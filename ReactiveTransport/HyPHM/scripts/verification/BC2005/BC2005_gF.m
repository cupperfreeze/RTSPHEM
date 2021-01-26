% Flux boundary condition acc. to
% C. Bahriawati, C. Carstensen
% ``Three Matlab Implementations of the Lowest-Order Raviart-Thomas
%  MFEM with a Posteriori Error Control''
% Copyright C. Bahriawati, C. Carstensen
%
% Input
%   t      [ 1 x 1 ] : (not provided)
%   x      [ 2 x 1 ] : evaluation point
%   specNo [ 1 x 1 ] : (not provided)
%   nu     [ 2 x 1 ] : outer unit normal
%
% Output
%   ret    [ 1 x 1 ] : flux boundary function gF
%
% Copyright 2009, 2010 F. Frank, Chair of Applied Mathematics I,
%           Department of Mathematics, University of Erlangen-Nuremberg,
%           G E R M A N Y

function ret = BC2005_gF(x, nu)

assert(isequal(size(x), [2, 1]))
assert(isequal(size(nu), [2, 1]))

[phi, r] = cart2pol(x(1), x(2));
phi = mod(phi, 2*pi);

ret = -(2 / 3 * r^(-1 / 3) * [-sin(phi / 3), cos(phi / 3)]) * nu;

return

end
