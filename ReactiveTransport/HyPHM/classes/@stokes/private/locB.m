%> @file classes/@Stokes/locB.m Performs the assembly of a local assembly matrix @f$\mathbf{B}_T@f$ under consideration of the global orientation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%>Performs the assembly of a local assembly matrix @f$\mathbf{B}_T@f$ under consideration of the global orientation.
%>
%> @param  ignoredArg   [Stokes].
%> @param  DF           Jacobian of affine mapping @f$F:\hat{T}\rightarrow T@f$ [2 x 2].
%> @retval ret          Local assembly matrix @f$\mathbf{B}_T@f$ [3 x 6].
%>
%>   @f[\mathbf{B}_{N'N} = (\psi_{N'}, \partial_x \phi_{N}).@f]
%>


function stiffB = locB(~, DF)

B1 = [-1 / 6, 0, 0, 1 / 6, -1 / 6, 1 / 6; ...
    0, 1 / 6, 0, 1 / 6, -1 / 6, -1 / 6; ...
    0, 0, 0, 1 / 3, -1 / 3, 0];

B2 = [-1 / 6, 0, 0, 1 / 6, 1 / 6, -1 / 6; ...
    0, 0, 0, 1 / 3, 0, -1 / 3; ...
    0, 0, 1 / 6, 1 / 6, -1 / 6, -1 / 6];

stiffB = DF(2, 2) * B1 - DF(2, 1) * B2;

end
