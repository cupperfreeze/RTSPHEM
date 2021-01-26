% Right hand side acc. to
% Florin A. Radu, Markus Bause, Alexander Prechtel and Sabine Attinger
% ``A mixed hybrid finite element discretization scheme for reactive transport
% in porous media''. In Numerical Mathematics and Advanced Applications, K.
% Kunisch, G. Of, O. Steinbach (editors), Springer, 2008, pp. 513-520.
%
% WITH ONE SPECIES ONLY AND NO REACTION TERM
%
% Input
%   t      [ 1     ] : evaluation time
%   x      [ 2 x 1 ] : evaluation point
%   specNo [ 1     ] : number of species, 1 for acceptor, 2 for donator
%
% Output
%   ret    [ 1     ] : coefficient function F
%

function ret = RBPA2008_F_onespec(t, x, specNo)

if ~isequal(size(x), [2, 1])
    error('Wrong arguments.')
end

y = x(2);
x = x(1);

if specNo == 1 % acceptor
    ret = -0.1111111111e-1 * (x - 1)^2 * y^2 * exp(-.1*t) ...
        -0.2222222222e-1 * y^2 * exp(-.1*t) ...
        -0.2222222222e-1 * (x - 1)^2 * exp(-.1*t) ...
        -2 / 9 * (x - 1)^2 * y * exp(-.1*t);
else
    error('Requested species not available.')
end

return

end
