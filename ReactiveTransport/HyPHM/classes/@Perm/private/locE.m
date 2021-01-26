%> @file classes/@Perm/locE.m Performs the assembly of a local assembly matrix @f$\mathbf{E}_T@f$ under consideration of the global orientation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%>Performs the assembly of a local assembly matrix @f$\mathbf{E}_T@f$ under consideration of the global orientation.
%>
%> @param  ignoredArg   [Perm].
%> @param  DF           Jacobian of affine mapping @f$F:\hat{T}\rightarrow T@f$ [numT x 2 x 2].
%> @retval ret          Local assembly matrix vs constant basis vector [numT x 1 x 6].
%>
%>   @f[\mathbf{E}_{N'N} = (\phi_N, \phi_{N'}).@f]
%>


function ret = locE(~, DF, numT)


% ret =  det(DF) * ...
%   [0, 0, 0, 1/6, 1/6, 1/6];
ret = (DF(:, 1, 1) .* DF(:, 2, 2) - DF(:, 1, 2) .* DF(:, 2, 1)) .* [0, 0, 0, 1 / 6, 1 / 6, 1 / 6];
end
