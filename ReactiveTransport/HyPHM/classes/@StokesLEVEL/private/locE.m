%> @file classes/@StokesLEVEL/locE.m Performs the assembly of a local assembly matrix @f$\mathbf{E}_T@f$ under consideration of the global orientation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%>Performs the assembly of a local assembly matrix @f$\mathbf{E}_T@f$ under consideration of the global orientation.
%>
%> @param  ignoredArg   [StokesLEVEL].
%> @param  DF           Jacobian of affine mapping @f$F:\hat{T}\rightarrow T@f$ [2 x 2].
%> @retval ret          Local assembly matrix @f$\mathbf{E}_T@f$ [6 x 6].
%>
%>   @f[\mathbf{E}_{N'N} = (\phi_N, \phi_{N'}).@f]
%>


function ret = locE(~, DF)

ret = det(DF) * ...
    [1 / 60, -1 / 360, -1 / 360, -1 / 90, 0, 0; ...
    -1 / 360, 1 / 60, -1 / 360, 0, -1 / 90, 0; ...
    -1 / 360, -1 / 360, 1 / 60, 0, 0, -1 / 90; ...
    -1 / 90, 0, 0, 4 / 45, 2 / 45, 2 / 45; ...
    0, -1 / 90, 0, 2 / 45, 4 / 45, 2 / 45; ...
    0, 0, -1 / 90, 2 / 45, 2 / 45, 4 / 45];

end
