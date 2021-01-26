%> @file +P1/setfh3data.m Set data from f = @(t, X, Y), where t [scalar], X [matrix], Y [matrix].
function ret = setfh3data(time, grid, fun3)

ret = fun3(time, grid.coordV(:, 1), grid.coordV(:, 2));

end
