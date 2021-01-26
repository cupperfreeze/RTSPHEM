%> @file classes/@StokesLEVEL/locA.m Performs the assembly of a local assembly matrix @f$\mathbf{A}_T@f$ under consideration of the global orientation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%>Performs the assembly of a local assembly matrix @f$\mathbf{A}_T@f$ under consideration of the global orientation.
%>
%> @param  ignoredArg   [StokesLEVEL].
%> @param  DF           Jacobian of affine mapping @f$F:\hat{T}\rightarrow T@f$ [numT x 2 x 2].
%> @retval ret          Local assembly matrix @f$\mathbf{A}_T@f$ [numT x 6 x 6].
%>
%>   @f[\mathbf{A}_{N'N} = (\vec{\nabla}\phi_N,\vec{\nabla}\phi_{N'}).@f]
%>


function ret = locA(~, DF, numT)


A1 = ... % dxdx
    [1 / 2, 1 / 6, 0, 0, 0, -2 / 3; ...
    1 / 6, 1 / 2, 0, 0, 0, -2 / 3; ...
    0, 0, 0, 0, 0, 0; ...
    0, 0, 0, 4 / 3, -4 / 3, 0; ...
    0, 0, 0, -4 / 3, 4 / 3, 0; ...
    -2 / 3, -2 / 3, 0, 0, 0, 4 / 3];


A2 = ... % dxdy + dydx
    [1, 1 / 6, 1 / 6, 0, -2 / 3, -2 / 3; ...
    1 / 6, 0, -1 / 6, 2 / 3, 0, -2 / 3; ...
    1 / 6, -1 / 6, 0, 2 / 3, -2 / 3, 0; ...
    0, 2 / 3, 2 / 3, 4 / 3, -4 / 3, -4 / 3; ...
    -2 / 3, 0, -2 / 3, -4 / 3, 4 / 3, 4 / 3; ...
    -2 / 3, -2 / 3, 0, -4 / 3, 4 / 3, 4 / 3];


A3 = ... % dydy
    [1 / 2, 0, 1 / 6, 0, -2 / 3, 0; ...
    0, 0, 0, 0, 0, 0; ...
    1 / 6, 0, 1 / 2, 0, -2 / 3, 0; ...
    0, 0, 0, 4 / 3, 0, -4 / 3; ...
    -2 / 3, 0, -2 / 3, 0, 4 / 3, 0; ...
    0, 0, 0, -4 / 3, 0, 4 / 3];


% G1 =  dot(DF(:, 2), DF(:, 2));
% G2 = -dot(DF(:, 1), DF(:, 2));
% G3 =  dot(DF(:, 1), DF(:, 1));
%ret = (G1*A1 + G2*A2 + G3*A3)/det(DF);
G1 = reshape(DF(:, 1, 2).*DF(:, 1, 2)+DF(:, 2, 2).*DF(:, 2, 2), [1, 1, numT]);
G2 = -reshape(DF(:, 1, 1).*DF(:, 1, 2)+DF(:, 2, 1).*DF(:, 2, 2), [1, 1, numT]);
G3 = reshape(DF(:, 1, 1).*DF(:, 1, 1)+DF(:, 2, 1).*DF(:, 2, 1), [1, 1, numT]);

k = reshape((DF(:, 1, 1).*DF(:, 2, 2)-DF(:, 2, 1).*DF(:, 1, 2)), 1, 1, numT);
ret = (G1 .* A1 + G2 .* A2 + G3 .* A3) ./ k;
end
