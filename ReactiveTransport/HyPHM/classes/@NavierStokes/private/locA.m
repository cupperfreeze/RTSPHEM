% locA performs the assembly of a local assembly matrix
% under consideration of the global orientation.
%
% Input
%   d        [ NavierStokes ]  data required by solver
%   kT       [ scalar ]        number of current element
%
% Output
%   stiff    [6 x 6]             local assembly matrix
%
%   stiff(N', N) = <phi_N, phi_N'> ,  N, N' = 1,...,6.
%
%
% See also NavierStokes.computeLevel
%
% Copyright 2009, 2010 F. Frank, Chair of Applied Mathematics I,
%           Department of Mathematics, University of Erlangen-Nuremberg,
%           G E R M A N Y

function stiff = locA(d, kT)


% analytical computation shows that (see @NavierStokes/algebra)
stiff = [1 / 60, -1 / 360, -1 / 360, -1 / 90, 0, 0; ...
    -1 / 360, 1 / 60, -1 / 360, 0, -1 / 90, 0; ...
    -1 / 360, -1 / 360, 1 / 60, 0, 0, -1 / 90; ...
    -1 / 90, 0, 0, 4 / 45, 2 / 45, 2 / 45; ...
    0, -1 / 90, 0, 2 / 45, 4 / 45, 2 / 45; ...
    0, 0, -1 / 90, 2 / 45, 2 / 45, 4 / 45];

stiff = 2 * d.grid.areaT(kT) * stiff;

end

% Remark: d.isEvolution is the flag which decides whether the evolution
% term should be taken into account.
