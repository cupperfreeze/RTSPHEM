% Visualization of the flux solution of BC2005.
%
% Copyright 2009 F. Frank, Chair of Applied Mathematics I,
%           Department of Mathematics, University of Erlangen-Nuremberg,
%           G E R M A N Y

function BC2005_q_disp()

[X, Y] = meshgrid(-1:.1:1);
n = length(X);
Z = zeros(n, n, 2);
for k = 1:n
    for j = 1:n
        Z(k, j, :) = BC2005_q(0, [X(k, j); Y(k, j)], 1);
    end
end

quiver(X, Y, Z(:, :, 1), Z(:, :, 2));

end
