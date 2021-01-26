% Copyright 2009, 2010, 2012 F. Frank, Chair of Applied Mathematics I,
%           Department of Mathematics, University of Erlangen-Nuremberg,
%           G E R M A N Y
%
% remark: All computations seem now to be speed-optimized by avoiding all
%         long for-loops.

function [coordV0T, baryE0T] = evalLocGeometry(coordV, V0T, E0T, V0E)

numT = size(V0T, 1);

% coordV0T (deps.: VOT)
coordV0T = zeros(numT, 3, 2);
for k = 1:3
    coordV0T(:, k, 1) = coordV(V0T(:, k), 1)';
    coordV0T(:, k, 2) = coordV(V0T(:, k), 2)';
end


% baryE0T (deps.: EOT, (baryE))
baryE = 0.5 * (coordV(V0E(:, 1), :) + coordV(V0E(:, 2), :));

baryE0T = zeros(numT, 3, 2);
for k = 1:3
    baryE0T(:, k, 1) = baryE(E0T(:, k), 1)';
    baryE0T(:, k, 2) = baryE(E0T(:, k), 2)';
end

end % evalLocGeometry
