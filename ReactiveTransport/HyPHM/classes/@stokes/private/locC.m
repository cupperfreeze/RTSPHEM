%> @file classes/@Stokes/locC.m Performs the assembly of a local assembly matrix @f$\mathbf{C}_T@f$ under consideration of the global orientation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%>Performs the assembly of a local assembly matrix @f$\mathbf{C}_T@f$ under consideration of the global orientation.
%>
%> @param  ignoredArg   [Stokes].
%> @param  DF           Jacobian of affine mapping @f$F:\hat{T}\rightarrow T@f$ [2 x 2].
%> @retval ret          Local assembly matrix @f$\mathbf{C}_T@f$ [3 x 6].
%>
%>   @f[\mathbf{C}_{N'N} = (\psi_{N'}, \partial_y \phi_{N}).@f]
%>

function stiff = locC(~, DF)

C1 = [-1 / 6, 0, 0, 1 / 6, -1 / 6, 1 / 6; ...
    0, 1 / 6, 0, 1 / 6, -1 / 6, -1 / 6; ...
    0, 0, 0, 1 / 3, -1 / 3, 0];

C2 = [-1 / 6, 0, 0, 1 / 6, 1 / 6, -1 / 6; ...
    0, 0, 0, 1 / 3, 0, -1 / 3; ...
    0, 0, 1 / 6, 1 / 6, -1 / 6, -1 / 6];

stiff = -DF(1, 2) * C1 + DF(1, 1) * C2;

end
