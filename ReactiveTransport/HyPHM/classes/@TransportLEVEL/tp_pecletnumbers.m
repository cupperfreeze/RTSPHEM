%> @file tp_pecletnumbers Calculation of the local Peclet numbers (aka grid Peclet numbers).
%>
%> The <i>local Peclet number</i> or <i>grid Peclet number</i> for a triangle @f$T\in\mathfrak{T}@f$ is defined by
%>
%> @f[\mathrm{Pe}_T = \frac{\|\vec{C}\|_{L^\infty(T)}\,h_T}{2\,\|\mathbf{D}\|_{L^\infty(T)} }.@f]
%>
%> This dimensionless number describes the ratio of the local advective to the local diffusive transport rate.
%> If @f$\mathrm{Pe}_T>1@f$ the convective part dominates the diffusive one and the flux is said to be (locally) <i>convection dominated</i> (and vice versa for @f$\mathrm{Pe}_T\leqq 1@f$).
%> Obviously, the local Peclet number can be decreased by using finer
%> grids.


function Pe = tp_pecletnumbers(d, dataC, dataD)

g = d.grid;
dataC_P0P0 = RT0.RT0toP0P0slice(g, dataC);
Pe = 0.5 * g.areaT .* max(abs(dataC_P0P0), [], 2) ./ max(abs(dataD), [], 2);

end
