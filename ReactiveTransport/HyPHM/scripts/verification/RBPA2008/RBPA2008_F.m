%> @file RBPA2008_F.m Right hand side acc. to @ref RBPA2008.

%
% Input
%   t      [ 1     ] : evaluation time
%   x      [ 2 x 1 ] : evaluation point
%   specNo [ 1     ] : number of species, 1 for acceptor, 2 for donator
%
% Output
%   ret    [ 1     ] : coefficient function F
%

function ret = RBPA2008_F(t, x, specNo)

if ~isequal(size(x), [2, 1])
    error('Wrong arguments.')
end

y = x(2);
x = x(1);

if specNo == 1 % acceptor
    ret = -0.1111111111e-1 * (x - 1)^2 * y^2 * exp(-.1*t) - 0.2222222222e-1 * y^2 * exp(-.1*t) ...
        -0.2222222222e-1 * (x - 1)^2 * exp(-.1*t) - 2 / 9 * (x - 1)^2 * y * exp(-.1*t) ...
        +2 / 2187 * (x - 1)^4 * y^7 * (exp(-.1 * t))^3 * x * (2 - x);

elseif specNo == 2 % donator
    ret = -0.3703703704e-2 * x * (2 - x) * y^3 * exp(-.1*t) + 0.7407407407e-2 * y^3 * exp(-.1*t) ...
        -0.2222222222e-1 * x * (2 - x) * y * exp(-.1*t) - 1 / 9 * x * (2 - x) * y^2 * exp(-.1*t) ...
        +1 / 2187 * (x - 1)^4 * y^7 * (exp(-.1 * t))^3 * x * (2 - x);
else
    error('Requested species not available.')
end

return

end
