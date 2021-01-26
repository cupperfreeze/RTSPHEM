% locF performs the assembly of a local assembly matrix
% under consideration of the global orientation.
%
% Input
%   d        [ NavierStokes ]  data required by solver
%   kT       [ scalar ]        number of current element
%
% Output
%   stiff    [3 x 6]             local assembly matrix
%
%   stiff(N', N) = <lambda_N', dxphi_N> ,   N' = 1,...,3, N = 1,...,6.
%
%
% See also NavierStokes.computeLevel
%
% Copyright 2009, 2010 F. Frank, Chair of Applied Mathematics I,
%           Department of Mathematics, University of Erlangen-Nuremberg,
%           G E R M A N Y

function stiff = locF(d, kT)

F1 = [-1 / 6, 0, 0, 1 / 6, -1 / 6, 1 / 6; ...
    0, 1 / 6, 0, 1 / 6, -1 / 6, -1 / 6; ...
    0, 0, 0, 1 / 3, -1 / 3, 0];

F2 = [-1 / 6, 0, 0, 1 / 6, 1 / 6, -1 / 6; ...
    0, 0, 0, 1 / 3, 0, -1 / 3; ...
    0, 0, 1 / 6, 1 / 6, -1 / 6, -1 / 6];

A = squeeze(d.grid.A(kT, :, :));

stiff = A(2, 2) * F1 - A(2, 1) * F2;

end

% Remark: d.isEvolution is the flag which decides whether the evolution
% term should be taken into account.