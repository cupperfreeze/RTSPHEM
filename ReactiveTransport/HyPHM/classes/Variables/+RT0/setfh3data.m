%> @file +RT0/setfh3data.m Set data from f = @(t, X, Y) {U(t,X,Y), V(t,X,Y)}, where t [scalar], X [matrix], Y [matrix].
%> Set data from <tt>f = @(t, X, Y) {U(t,X,Y), V(t,X,Y)}</tt>, where t [scalar], X
%> [matrix], Y [matrix] and U, V the components of f in cartesian coordinates which can take
%> matrices X and Y.  Note that the target of f must be of type
%> <tt>cell</tt>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ret = setfh3data(time, grid, fun3)

UV_cart = fun3(time, grid.baryE(:, 1), grid.baryE(:, 2));
ret = UV_cart{1} .* grid.nuE(:, 1) + UV_cart{2} .* grid.nuE(:, 2);

end
