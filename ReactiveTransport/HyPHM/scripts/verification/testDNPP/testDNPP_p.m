%> @file testDNPP_p.m Exact pressure distribution @f$p@f$
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%>
%> @param t evaluation time [scalar]
%> @param x sample point in domain @f$[2\times 1]@f$
%> @retval ret see function description [scalar]

function ret = testDNPP_p(t, x)
msg = 'HyPHM: RTFM.';

assert(isequal(size(x), [2, 1]), msg)
assert(isscalar(t), msg)

ret = -1 / 4 * (cos(2 * pi * x(1)) + cos(2 * pi * x(2)));

end
