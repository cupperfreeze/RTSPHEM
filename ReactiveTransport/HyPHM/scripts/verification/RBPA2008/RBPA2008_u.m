% Analytical solution of concentration acc. to
% Florin A. Radu, Markus Bause, Alexander Prechtel and Sabine Attinger
% ``A mixed hybrid finite element discretization scheme for reactive transport
% in porous media''. In Numerical Mathematics and Advanced Applications, K.
% Kunisch, G. Of, O. Steinbach (editors), Springer, 2008, pp. 513-520.
%
% Input
%   x        [ 2 x 1  ] : evaluation point
%   specNo   [ 1      ] : number of species, 1 for acceptor, 2 for donator
%
% Output
%   ret      [ 1      ] : c_i(t, x)
%
% Copyright 2009 F. Frank, Chair of Applied Mathematics I,
%           Department of Mathematics, University of Erlangen-Nuremberg,
%           G E R M A N Y

function ret = RBPA2008_u(t, x, specNo)

if ~isequal(size(x), [2, 1])
    error('Wrong arguments used in RBPA2008_u().')
end

y = x(2);
x = x(1);

if specNo == 1 % acceptor
    ret = 1 / 9 * (x - 1)^2 * y^2 * exp(-.1*t);
elseif specNo == 2 % donator
    ret = 1 / 27 * x * (2 - x) * y^3 * exp(-.1*t);
else
    error('Requested species not available.')
end

return

end
