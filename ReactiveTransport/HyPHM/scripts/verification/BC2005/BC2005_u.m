% Exact scalar solution acc. to
% C. Bahriawati, C. Carstensen
% ``Three Matlab Implementations of the Lowest-Order Raviart-Thomas
%  MFEM with a Posteriori Error Control''
%
% Input
%   t      [ 1 x 1 ] : (not provided)
%   x      [ 2 x 1 ] : evaluation point
%   specNo [ 1 x 1 ] : (not provided)
%
% Output
%   ret    [ 1 x 1 ] : u(t, x)
%
% Copyright 2009, 2010 F. Frank, Chair of Applied Mathematics I,
%           Department of Mathematics, University of Erlangen-Nuremberg,
%           G E R M A N Y

function ret = BC2005_u(~, x, specNo)

assert(isequal(size(x), [2, 1]), 'HyPHM: Wrong arguments used in BC2005_q().')
assert(specNo == 1, 'HyPHM: Requested species not available.')

[phi, r] = cart2pol(x(1), x(2));
phi = mod(phi, 2*pi);

ret = r^(2 / 3) * sin(2*phi/3);

return

end
